// # Literate Flatbush
// ## Understanding a fast, elegant RTree implementation
//
// [Kyle Barron](https://kylebarron.dev/) \
// January 1, 2025
//
// Spatial indexes are a tried and true part of geospatial software engineering that speed up queries. But ever wondered how an RTree is actually implemented?
//
// TODO: fix this intro.
// TODO: somewhere in intro explain what a hilbert curve is.
//
// This post is a ["literate"](https://en.wikipedia.org/wiki/Literate_programming) fork of Flatbush, a blazing-fast, memory-efficient RTree. Documentation and code are interspersed, letting you follow along with the code.
//
// The code in this particular post is in JavaScript, because that's what the canonical implementation uses. But it's the _algorithm_ that's important here, so don't get too caught up in the JavaScript.
//
// ## Overview
//
// The Flatbush algorithm generates a **static, packed, ABI-stable RTree**.
// Let's break that down:
//
// - [**RTree**](https://en.wikipedia.org/wiki/R-tree): a [spatial
//   index](https://blog.mapbox.com/a-dive-into-spatial-search-algorithms-ebd0c5e39d2a)
//   for storing geospatial vector data that allows for fast spatial queries.
//
//   It's a form of a
//   ["tree"](https://en.wikipedia.org/wiki/Tree_(abstract_data_type)), where
//   there's one root node that has `nodeSize` children. Each of those nodes
//   have their own `nodeSize` children, and so on. The tree structure allows
//   you to avoid superfluous checks and quickly find matching candidates for
//   your query. In particular, an RTree stores a _bounding box_ for each
//   geometry.
// - **static**: the index is immutable. All geometries need to be added to the index before any searches can be done. Geometries can't be added to an existing index later.
// - **packed**: all nodes are at full capacity (except for the last node at each tree level). Because the tree is static, we don't need to reserve space in each node for future additions. This improves memory efficiency.
// - **ABI-stable**: the underlying binary memory layout is stable. This enables zero-copy sharing between threads or, in my Rust port, between Rust and Python.
//
// There are a few really nice features about this algorithm, and why it makes
// sense to document it like this:
//
// - **Speed**: This is likely the fastest static spatial index in JavaScript. Ports of the algorithm are among the fastest spatial indexes in other languages, too.
// - **Single, contiguous underlying buffer**: The index is contained in a single `ArrayBuffer`, which makes it easy to share across multiple threads or persist and use later.
// - **Memory-efficiency**: because the index is fully packed, it's highly memory efficient.
// - **Bounded-memory**: for any given number of items and node size, you can infer the total memory that will be used by the RTree.
// - **Simple and elegant**: The JavaScript implementation is just a few hundred lines of code, and in my opinion it's quite elegant how the structure of the tree implicitly maintains the insertion index.
//
// Keep in mind there are a few restrictions as well:
//
// - Only two-dimensional data. Because the algorithm uses powers of two, only two-dimensional data is supported.
//
//
//
// I've written _two_ ports of Flatbush. I originally explored a (now-archived)
// [Cython port](https://github.com/kylebarron/pyflatbush) but deprecated that
// in favor of writing a [stable, fast, implementation in
// Rust](https://github.com/kylebarron/geo-index), with [Python
// bindings](https://github.com/kylebarron/geo-index/tree/main/python).
//
// This is a "literate" fork of the [`flatbush`][] library. Markdown-formatted
// comments are written in the code, and then
// [`docco`](https://ashkenas.com/docco/) converts the file into a rendered HTML
// file. No code modifications are made in this fork. Only comments are added.
//
// All credit for this code goes to Volodymyr Agafonkin and contributors to the
// [`flatbush`][] project, forked here under the ISC license. Any errors in
// explanation are mine alone.
//
// [`flatbush`]: https://github.com/mourner/flatbush

import FlatQueue from "flatqueue";

// Flatbush supports a variety of
// [`TypedArray`](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/TypedArray)
// types store box coordinate data. Flatbush uses `Float64Array` by default.
const ARRAY_TYPES = [
  Int8Array,
  Uint8Array,
  Uint8ClampedArray,
  Int16Array,
  Uint16Array,
  Int32Array,
  Uint32Array,
  Float32Array,
  Float64Array,
];

// The Flatbush serialized format version is bumped whenever the binary layout
// of the index changes
const VERSION = 3; // serialized format version

/** @typedef {Int8ArrayConstructor | Uint8ArrayConstructor | Uint8ClampedArrayConstructor | Int16ArrayConstructor | Uint16ArrayConstructor | Int32ArrayConstructor | Uint32ArrayConstructor | Float32ArrayConstructor | Float64ArrayConstructor} TypedArrayConstructor */

// ## Flatbush
//
// The `Flatbush` class is the only export from the Flatbush library. It
// contains functions to create and query the spatial index.
export default class Flatbush {
  // ### Flatbush.from
  //
  // One of Flatbush's goals is to support zero-copy usage, meaning that you can
  // take an `ArrayBuffer` backing a Flatbush index and
  // [_transfer_](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API/Transferable_objects)
  // it between threads much faster than copying it via a [structured
  // clone](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API/Structured_clone_algorithm).
  // (Transfering is `O(1)` while cloning is `O(n)`.)
  //
  // The `from` static method on the class reconstructs a `Flatbush` instance
  // from a raw `ArrayBuffer`.

  /**
   * Recreate a Flatbush index from raw `ArrayBuffer` or `SharedArrayBuffer` data.
   * @param {ArrayBuffer | SharedArrayBuffer} data
   * @returns {Flatbush} index
   */
  static from(data) {
    if (!data || data.byteLength === undefined || data.buffer) {
      throw new Error(
        "Data must be an instance of ArrayBuffer or SharedArrayBuffer."
      );
    }

    // The first 8 bytes of the Flatbush buffer contain a header with metadata
    // describing the contained index.
    //
    // - byte 1: a "magic byte" set to `0xfb`.
    // - byte 2: 4 bits for the serialized format version and four bits for the array type used for storing coordinates
    // - byte 3-4: a uint16-encoded number representing the size of each node
    // - byte 5-8: a uint32-encoded number representing the total number of
    //   items in the index.
    const [magic, versionAndType] = new Uint8Array(data, 0, 2);
    if (magic !== 0xfb) {
      throw new Error("Data does not appear to be in a Flatbush format.");
    }
    const version = versionAndType >> 4;
    if (version !== VERSION) {
      throw new Error(`Got v${version} data when expected v${VERSION}.`);
    }
    const ArrayType = ARRAY_TYPES[versionAndType & 0x0f];
    if (!ArrayType) {
      throw new Error("Unrecognized array type.");
    }
    const [nodeSize] = new Uint16Array(data, 2, 1);
    const [numItems] = new Uint32Array(data, 4, 1);

    // Given the above parsed metadata about the index, pass these to the
    // constructor. Because the `data` argument (passed last) is not
    // `undefined`, the constructor will not create a new underlying buffer, but
    // rather reuse the existing buffer.
    return new Flatbush(numItems, nodeSize, ArrayType, undefined, data);
  }

  // ### Constructor
  //
  // The Flatbush constructor initializes the memory space (`ArrayBuffer`) for a
  // Flatbush tree given the number of items the tree will contain and the
  // number of elements per tree node.
  /**
   * Create a Flatbush index that will hold a given number of items.
   * @param {number} numItems
   * @param {number} [nodeSize=16] Size of the tree node (16 by default).
   * @param {TypedArrayConstructor} [ArrayType=Float64Array] The array type used for coordinates storage (`Float64Array` by default).
   * @param {ArrayBufferConstructor | SharedArrayBufferConstructor} [ArrayBufferType=ArrayBuffer] The array buffer type used to store data (`ArrayBuffer` by default).
   * @param {ArrayBuffer | SharedArrayBuffer} [data] (Only used internally)
   */
  constructor(
    numItems,
    nodeSize = 16,
    ArrayType = Float64Array,
    ArrayBufferType = ArrayBuffer,
    data
  ) {
    if (numItems === undefined)
      throw new Error("Missing required argument: numItems.");
    if (isNaN(numItems) || numItems <= 0)
      throw new Error(`Unexpected numItems value: ${numItems}.`);

    this.numItems = +numItems;
    this.nodeSize = Math.min(Math.max(+nodeSize, 2), 65535);

    // This do-while loop calculates the total number of nodes at each level of
    // the R-tree (and thus also the total number of nodes). This will be used
    // to know the amount of bytes to allocate for each level of the tree.
    //
    // The tree is laid out in memory from bottom to top. `_levelBounds` is an
    // array that stores the offset within the coordinates array where each
    // level **ends**. The first element of `_levelBounds` is `n * 4`, meaning
    // that the slice of the coordinates array from `0` to `n * 4` contains the
    // bottom (leaves) of the tree.
    //
    // Then the slice of the coordinates array from `_levelBounds[0]` to
    // `_levelBounds[1]` represents the first level of the tree, that is, the
    // direct parent nodes of the leaves. And so on, `_levelBounds[1]` to
    // `_levelBounds[2]` represents the nodes at level 2, the grandparent nodes
    // of the leaf nodes.
    //
    // So for example if `numItems` is 10,000 and `nodeSize` is 16,
    // `levelBounds` will be:
    // ```
    // [40000, 42500, 42660, 42672, 42676]
    // ```
    //
    // That is:
    // - The first 40,000 elements are coordinates of the leaf nodes (4 coordinates per node). So there are 10,000 leaves.
    // - 2,500 coordinates and 625 nodes one level higher
    // - 160 coordinates and 40 nodes two levels higher
    // - 12 coordinates and 3 nodes three levels higher
    // - 1 root node four levels higher, at the top of the tree, with a single 4-coordinate box.
    //
    // Keep in mind that because this is a _packed_ tree, every node within a
    // single level will be **completely full** (contain exactly `nodeSize`
    // elements) except for the last node.
    //
    // `numNodes` ends up as the total number of nodes in the tree, including
    // all leaves.
    let n = numItems;
    let numNodes = n;
    this._levelBounds = [n * 4];
    do {
      n = Math.ceil(n / this.nodeSize);
      numNodes += n;
      this._levelBounds.push(numNodes * 4);
    } while (n !== 1);

    // Flatbush doesn't manage references to objects directly. Rather, it
    // operates in terms of the _insertion index_. So Flatbush only internally
    // maintains these insertion indices.
    //
    // `IndexArrayType` will be used to create the `indices` array, to store the
    // ordering of the input boxes. If possible, a `Uint16Array` will be used to
    // save space. If the values would overflow a `Uint16Array`, a `Uint32Array`
    // is used. The largest number a `Uint16Array` can hold is `2^16 = 65,536`.
    // Since each node holds four values, this gets divided by `4` and `65,536 /
    // 4 = 16,384`. This is why the check here is for 16,384.
    this.ArrayType = ArrayType;
    this.IndexArrayType = numNodes < 16384 ? Uint16Array : Uint32Array;

    // In order to accurately interpret the index from raw bytes, we need to
    // record in the header which index type we're using.
    const arrayTypeIndex = ARRAY_TYPES.indexOf(this.ArrayType);

    // The number of bytes needed to store all box coordinate data for all
    // nodes.
    const nodesByteSize = numNodes * 4 * this.ArrayType.BYTES_PER_ELEMENT;

    if (arrayTypeIndex < 0) {
      throw new Error(`Unexpected typed array class: ${ArrayType}.`);
    }

    // This `if` statement switches on whether the `data` argument was passed in
    // (i.e. this constructor is called by `Flatbush.from`). If `data` exists,
    // this will create the `_boxes` and `_indices` arrays as **views** on the
    // existing `ArrayBuffer` without allocating any new memory.
    if (data && data.byteLength !== undefined && !data.buffer) {
      this.data = data;

      // `this._boxes` is created as a view on `this.data` starting after the
      // header (8 bytes) and with `numNodes * 4` elements. `this._indices` is
      // created as a view on `this.data` starting after the end of
      // `this._boxes` and containing `numNodes` elements.
      this._boxes = new this.ArrayType(this.data, 8, numNodes * 4);
      this._indices = new this.IndexArrayType(
        this.data,
        8 + nodesByteSize,
        numNodes
      );

      // The coordinate data in the `_boxes` array is stored from the leaves up.
      // So the last box is the single node that contains all data. The index of
      // the last box is the four values in `_boxes` up to `numNodes * 4`.
      //
      // This sets the total bounds on the `Flatbush` instance to the extent of
      // that box.
      //
      // We also set `this._pos` as the total number of coordinates. `this._pos`
      // is a pointer into the `this._boxes` array, used while adding new boxes
      // to the instance. This also allows for inferring whether the `Flatbush`
      // instance has been "finished" (sorted) or not.
      //
      // If the instance has already been sorted, it's impermissible to add more
      // data. If the instance has not yet been sorted, query methods may not
      // be called.
      this._pos = numNodes * 4;
      this.minX = this._boxes[this._pos - 4];
      this.minY = this._boxes[this._pos - 3];
      this.maxX = this._boxes[this._pos - 2];
      this.maxY = this._boxes[this._pos - 1];

      // In the `else` case, a `data` buffer was not provided, so we need to
      // allocate data for the backing buffer.
      //
      // `this.data` is a new `ArrayBuffer` with space for the header plus all
      // box data plus all index data. Then `this._boxes` is created as a view
      // on `this.data` starting after the header and with `numNodes * 4`
      // elements. `this._indices` is created as a view on `this.data` starting
      // after the end of `this._boxes`.
    } else {
      this.data = new ArrayBufferType(
        8 + nodesByteSize + numNodes * this.IndexArrayType.BYTES_PER_ELEMENT
      );
      this._boxes = new this.ArrayType(this.data, 8, numNodes * 4);
      this._indices = new this.IndexArrayType(
        this.data,
        8 + nodesByteSize,
        numNodes
      );

      // We set `this._pos` to 0. This means that no boxes have
      // yet been added to the index, and it tells any query methods to throw
      // until `finish` has been called.
      this._pos = 0;

      // The RTree needs to maintain its total bounds (the global bounding box
      // of all values) in order to set the bounds for the hilbert space.
      //
      // We initialize these bounds to `Infinity` values that will be corrected
      // when adding data. The minimum x/y of any box will be less than positive
      // infinity and the maximum x/y of any box will be greater than negative
      // infinity. The `add()` call will adjust these bounds if necessary.
      this.minX = Infinity;
      this.minY = Infinity;
      this.maxX = -Infinity;
      this.maxY = -Infinity;

      // These lines set the header values with metadata from the instance.
      //
      // The first byte, `0xfb` is a "magic byte", used as basic validation that
      // this buffer is indeed a Flatbush index.
      //
      // Since `arrayTypeIndex` is known to have only 9 values, it doesn't need
      // to take up a a full byte. Here it shares a single byte with the
      // Flatbush format version.
      new Uint8Array(this.data, 0, 2).set([
        0xfb,
        (VERSION << 4) + arrayTypeIndex,
      ]);
      new Uint16Array(this.data, 2, 1)[0] = nodeSize;
      new Uint32Array(this.data, 4, 1)[0] = numItems;
    }

    // We initialize a [priority
    // queue](https://en.wikipedia.org/wiki/Priority_queue) used for
    // k-nearest-neighbors queries in the `neighbors` method.
    /** @type FlatQueue<number> */
    this._queue = new FlatQueue();
  }

  // ### Flatbush.Add
  /**
   * Add a given rectangle to the index.
   * @param {number} minX
   * @param {number} minY
   * @param {number} maxX
   * @param {number} maxY
   * @returns {number} A zero-based, incremental number that represents the newly added rectangle.
   */
  add(minX, minY, maxX, maxY) {
    // We need to know the insertion index of the box presently being added.
    //
    // In the constructor, `this._pos` is initialized to `0` and in each call
    // to `add()`, `this._pos` is incremented by `4`. Dividing `this._pos` by
    // `4` retrieves the 0-based index of the box about to be inserted.
    //
    // This bit shift:
    // ```js
    // this._pos >> 2
    // ```
    // is equivalent to
    // ```js
    // this._pos / 4
    // ```
    // but the bit shift is faster because it informs the JS engine that we
    // expect the output to be an integer.
    //
    // Because there are 4 values for each item, using `_pos` is an easy way to
    // infer the insertion index without having to maintain a separate counter.
    const index = this._pos >> 2;
    const boxes = this._boxes;

    // We set the value of `this._indices` at the current index's position to
    // the value of the current index. So `this._indices` stores the insertion
    // index of each box.
    //
    // Later, inside the `finish` method, we'll sort the boxes by their hilbert
    // value and jointly reorder the values in `_indices`, ensuring that we keep
    // the indices and boxes in sync.
    //
    // This means that for any box representing a leaf node at position `i`
    // (where `i` points to a _box_ not a _coordinate_ inside a box),
    // `this._indices[i]` retrieves the original insertion-order index of that
    // box.
    this._indices[index] = index;

    // We set the coordinates of this box into the `boxes` array. Note that
    // `this._pos++` is evaluated **after** the box index is set. So
    //
    // ```js
    // boxes[this._pos++] = minX;
    // ```
    //
    // is equivalent to
    //
    // ```js
    // boxes[this._pos] = minX;
    // this._pos += 1;
    // ```
    boxes[this._pos++] = minX;
    boxes[this._pos++] = minY;
    boxes[this._pos++] = maxX;
    boxes[this._pos++] = maxY;

    // Update the total bounds of this instance if this rectangle is larger than
    // the existing bounds.
    if (minX < this.minX) this.minX = minX;
    if (minY < this.minY) this.minY = minY;
    if (maxX > this.maxX) this.maxX = maxX;
    if (maxY > this.maxY) this.maxY = maxY;

    return index;
  }

  // ### Flatbush.finish
  //
  // TODO: explain high level of sorting by hilbert curve.
  //
  /** Perform indexing of the added rectangles. */
  finish() {
    // Recall that in the `add` method, we increment `this._pos` by `1` for each
    // coordinate of each box. Here we validate that we've added the same number
    // of boxes as we provisioned in the constructor. Remember that `>> 2` is
    // equivalent to `/ 4`.
    if (this._pos >> 2 !== this.numItems) {
      throw new Error(
        `Added ${this._pos >> 2} items when expected ${this.numItems}.`
      );
    }
    const boxes = this._boxes;

    // If the total number of items in the tree is less than the node size, that
    // means we'll only have a single non-leaf node in the tree. In that case,
    // we don't even need to sort by hilbert value. We can just assign the total
    // bounds of the tree to the following box and return.
    if (this.numItems <= this.nodeSize) {
      boxes[this._pos++] = this.minX;
      boxes[this._pos++] = this.minY;
      boxes[this._pos++] = this.maxX;
      boxes[this._pos++] = this.maxY;
      return;
    }

    // Using the total bounds of the tree, we compute the height and width of
    // the hilbert space and instantiate space for the hilbert values.
    const width = this.maxX - this.minX || 1;
    const height = this.maxY - this.minY || 1;
    const hilbertValues = new Uint32Array(this.numItems);
    const hilbertMax = (1 << 16) - 1;

    // Map box centers into Hilbert coordinate space and calculate Hilbert
    // values using the `hilbert` function defined below.
    //
    // This for loop iterates over every box. At the beginning of each loop
    // iteration, `pos` is equal to `i * 4`.
    for (let i = 0, pos = 0; i < this.numItems; i++) {
      const minX = boxes[pos++];
      const minY = boxes[pos++];
      const maxX = boxes[pos++];
      const maxY = boxes[pos++];
      const x = Math.floor(
        (hilbertMax * ((minX + maxX) / 2 - this.minX)) / width
      );
      const y = Math.floor(
        (hilbertMax * ((minY + maxY) / 2 - this.minY)) / height
      );
      hilbertValues[i] = hilbert(x, y);
    }

    // Up until this point, the values in `boxes` and in `this._indices` are
    // still in _insertion order_. We need to jointly sort the boxes and indices
    // according to their hilbert values, which we do in this `sort` utility,
    // documented below.
    //
    // Note that this canonical implementation chooses to sort based on hilbert
    // value, but that's actually not necessary to maintain ABI-stability. Any
    // two-dimensional sort will work. My [Rust
    // port](https://github.com/kylebarron/geo-index) defines an [extensible
    // trait](https://docs.rs/geo-index/latest/geo_index/rtree/sort/trait.Sort.html)
    // for sorting and provides an implementation of [sort-tile-recursive
    // sorting](https://ia600900.us.archive.org/27/items/nasa_techdoc_19970016975/19970016975.pdf).
    sort(
      hilbertValues,
      boxes,
      this._indices,
      0,
      this.numItems - 1,
      this.nodeSize
    );

    // Now the leaves of the tree have been sorted, but we still need to
    // construct the rest of the tree.
    //
    // For each level of the tree, we need to generate parent nodes that contain
    // `nodeSize` child nodes. We do this starting from the leaves, working from
    // the bottom up.
    //
    // TODO: resume careful proofreading from here.
    //
    // Here the iteration variable, `i`, refers to the index into the
    // `this._levelBounds` array, which is also the positional level of the
    // tree. So when `i` is `0`, we're iterating over the original geometry
    // boxes. When `i` is `1`, we're iterating over the parent nodes one level
    // up that we previously generated from the first loop iteration.
    //
    // As elsewhere, `pos` is a local variable that points to a coordinate
    // within a box at the given level `i` of the tree.
    for (let i = 0, pos = 0; i < this._levelBounds.length - 1; i++) {
      // Next, we want to scan through all nodes at this level of the tree,
      // generating a parent node for each block of consecutive `nodeSize`
      // nodes.
      //
      // Here, `end` is index of the first coordinate at the _next level above
      // the current level_. Therefore, all index values in the range up until
      // but excluding `end` refer to box coordinates at the current level of
      // the tree.
      //
      // We then scan over all of these box coordinates in this while loop.
      const end = this._levelBounds[i];
      while (pos < end) {
        // The node index (which will later get set in the `indices` array) is
        // defined to be the pointer to the first element of a box.
        const nodeIndex = pos;

        // Calculate the bounding box for the new parent node.
        //
        // We initialize the bounding box to the first box and then expand the
        // box while looping over the rest of the elements that together are the
        // children of this parent node we're creating.
        //
        // Note the `j = 1` in the loop; this is a small optimization because we
        // initialize the `node*` variables to the first element.
        //
        // Also note that in the loop we constrain the iteration variable `j` to
        // be both less than the node size and for `pos < end`. The former
        // ensures we have only a maximum of `nodeSize` elements informing the
        // parent node's boundary. The latter ensures that we don't accidentally
        // overflow the current tree level.
        let nodeMinX = boxes[pos++];
        let nodeMinY = boxes[pos++];
        let nodeMaxX = boxes[pos++];
        let nodeMaxY = boxes[pos++];
        for (let j = 1; j < this.nodeSize && pos < end; j++) {
          nodeMinX = Math.min(nodeMinX, boxes[pos++]);
          nodeMinY = Math.min(nodeMinY, boxes[pos++]);
          nodeMaxX = Math.max(nodeMaxX, boxes[pos++]);
          nodeMaxY = Math.max(nodeMaxY, boxes[pos++]);
        }

        // Now that we know the extent of the parent node, we can add the new
        // node's information to the tree data.
        //
        // **TODO**: think a bit more/describe how the
        //
        // ```js
        // this._indices[this._pos >> 2] = nodeIndex;
        // ```
        //
        // line works... in particular how is `_indices` laid out across levels?
        // Can you create a diagram of this?
        //
        // Note that we're setting the parent node into `this._indices` and
        // `boxes` according to `this._pos`, which **is a different variable
        // than the local `pos` variable that's incremented in this loop.**
        //
        // Impressively, these loops do all the hard work of constructing the
        // tree! That's it! The structure of the tree and the coordinates of all
        // the parent nodes are now fully contained within `this._indices` and
        // `boxes`, which are both views on `this.data`!
        this._indices[this._pos >> 2] = nodeIndex;
        boxes[this._pos++] = nodeMinX;
        boxes[this._pos++] = nodeMinY;
        boxes[this._pos++] = nodeMaxX;
        boxes[this._pos++] = nodeMaxY;
      }
    }
  }

  // ### Flatbush.search
  /**
   * Search the index by a bounding box.
   * @param {number} minX
   * @param {number} minY
   * @param {number} maxX
   * @param {number} maxY
   * @param {(index: number) => boolean} [filterFn] An optional function for filtering the results.
   * @returns {number[]} An array of indices of items intersecting or touching the given bounding box.
   */
  search(minX, minY, maxX, maxY, filterFn) {
    // A simple check to ensure that this index has been finished/sorted.
    if (this._pos !== this._boxes.length) {
      throw new Error("Data not yet indexed - call index.finish().");
    }

    // `nodeIndex` is initialized to the root node. The root node is the last
    // node in `this._boxes`. We subtract `4` so that `nodeIndex` points to the
    // first coordinate of the box.
    //
    // `queue` holds integers that represent the `pos` of internal nodes that
    // still need to be searched. That is, where the parent node of the node
    // represented by `pos` intersected the search predicate.
    //
    // `results` holds integers that represent the insertion indexes of index rows
    // that match the search predicate.
    /** @type number | undefined */
    let nodeIndex = this._boxes.length - 4;
    const queue = [];
    const results = [];

    // Now we have our search loop.
    //
    // `while (nodeIndex !== undefined)` will be `true` as long as there are
    // still elements remaining in `queue` (note that the last line of the
    // `while` loop is `nodeIndex = queue.pop();`).
    while (nodeIndex !== undefined) {
      // find the end index of the node
      const end = Math.min(
        nodeIndex + this.nodeSize * 4,
        upperBound(nodeIndex, this._levelBounds)
      );

      // search through child nodes
      for (let /** @type number */ pos = nodeIndex; pos < end; pos += 4) {
        // check if node bbox intersects with query bbox
        if (maxX < this._boxes[pos]) continue; // maxX < nodeMinX
        if (maxY < this._boxes[pos + 1]) continue; // maxY < nodeMinY
        if (minX > this._boxes[pos + 2]) continue; // minX > nodeMaxX
        if (minY > this._boxes[pos + 3]) continue; // minY > nodeMaxY

        // `pos` is a pointer to the first coordinate of a given node.
        // `pos` has multiple purposes. In the
        //
        // `pos >> 2` is a faster way of expressing `pos / 4`, where we can
        // inform the JS engine that the output will be an integer.
        //
        // - If the current node _is not_ a leaf node, `index` is the `pos` of
        //   the first child node. This is a node that we should evaluate later,
        //   so we add it to the `queue` array.
        // - If the current node _is_ a leaf node, then `index` is the original
        //   insertion index, and we add it to the `results` array.
        //
        // Each box has 4 coordinates, and recall that the
        //
        // I believe `| 0` is just a JS engine optimization.
        const index = this._indices[pos >> 2] | 0;

        if (nodeIndex >= this.numItems * 4) {
          queue.push(index); // node; add it to the search queue
        } else if (filterFn === undefined || filterFn(index)) {
          results.push(index); // leaf item
        }
      }

      nodeIndex = queue.pop();
    }

    return results;
  }

  /**
   * Search items in order of distance from the given point.
   * @param {number} x
   * @param {number} y
   * @param {number} [maxResults=Infinity]
   * @param {number} [maxDistance=Infinity]
   * @param {(index: number) => boolean} [filterFn] An optional function for filtering the results.
   * @returns {number[]} An array of indices of items found.
   */
  neighbors(x, y, maxResults = Infinity, maxDistance = Infinity, filterFn) {
    // A simple check to ensure that this index has been finished/sorted.
    if (this._pos !== this._boxes.length) {
      throw new Error("Data not yet indexed - call index.finish().");
    }

    /** @type number | undefined */
    let nodeIndex = this._boxes.length - 4;
    const q = this._queue;
    const results = [];
    const maxDistSquared = maxDistance * maxDistance;

    outer: while (nodeIndex !== undefined) {
      // find the end index of the node
      const end = Math.min(
        nodeIndex + this.nodeSize * 4,
        upperBound(nodeIndex, this._levelBounds)
      );

      // add child nodes to the queue
      for (let pos = nodeIndex; pos < end; pos += 4) {
        const index = this._indices[pos >> 2] | 0;

        const dx = axisDist(x, this._boxes[pos], this._boxes[pos + 2]);
        const dy = axisDist(y, this._boxes[pos + 1], this._boxes[pos + 3]);
        const dist = dx * dx + dy * dy;
        if (dist > maxDistSquared) continue;

        if (nodeIndex >= this.numItems * 4) {
          q.push(index << 1, dist); // node (use even id)
        } else if (filterFn === undefined || filterFn(index)) {
          q.push((index << 1) + 1, dist); // leaf item (use odd id)
        }
      }

      // pop items from the queue
      while (q.length && q.peek() & 1) {
        const dist = q.peekValue();
        if (dist > maxDistSquared) break outer;
        results.push(q.pop() >> 1);
        if (results.length === maxResults) break outer;
      }

      nodeIndex = q.length ? q.pop() >> 1 : undefined;
    }

    q.clear();
    return results;
  }
}

// The remaining code is "just" utility functions.
//
// I won't document these because they tend to be self explanatory and this post
// is focused more on the RTree implementation itself.
/**
 * 1D distance from a value to a range.
 * @param {number} k
 * @param {number} min
 * @param {number} max
 */
function axisDist(k, min, max) {
  return k < min ? min - k : k <= max ? 0 : k - max;
}

/**
 * Binary search for the first value in the array bigger than the given.
 * @param {number} value
 * @param {number[]} arr
 */
function upperBound(value, arr) {
  let i = 0;
  let j = arr.length - 1;
  while (i < j) {
    const m = (i + j) >> 1;
    if (arr[m] > value) {
      j = m;
    } else {
      i = m + 1;
    }
  }
  return arr[i];
}

// `sort`: Custom quicksort that partially sorts bbox data alongside the hilbert values.
/**
 * Custom quicksort that partially sorts bbox data alongside the hilbert values.
 * @param {Uint32Array} values
 * @param {InstanceType<TypedArrayConstructor>} boxes
 * @param {Uint16Array | Uint32Array} indices
 * @param {number} left
 * @param {number} right
 * @param {number} nodeSize
 */
function sort(values, boxes, indices, left, right, nodeSize) {
  if (Math.floor(left / nodeSize) >= Math.floor(right / nodeSize)) return;

  const pivot = values[(left + right) >> 1];
  let i = left - 1;
  let j = right + 1;

  while (true) {
    do i++;
    while (values[i] < pivot);
    do j--;
    while (values[j] > pivot);
    if (i >= j) break;
    swap(values, boxes, indices, i, j);
  }

  sort(values, boxes, indices, left, j, nodeSize);
  sort(values, boxes, indices, j + 1, right, nodeSize);
}

// `swap`: Swap two values and two corresponding boxes.
/**
 * Swap two values and two corresponding boxes.
 * @param {Uint32Array} values
 * @param {InstanceType<TypedArrayConstructor>} boxes
 * @param {Uint16Array | Uint32Array} indices
 * @param {number} i
 * @param {number} j
 */
function swap(values, boxes, indices, i, j) {
  const temp = values[i];
  values[i] = values[j];
  values[j] = temp;

  const k = 4 * i;
  const m = 4 * j;

  const a = boxes[k];
  const b = boxes[k + 1];
  const c = boxes[k + 2];
  const d = boxes[k + 3];
  boxes[k] = boxes[m];
  boxes[k + 1] = boxes[m + 1];
  boxes[k + 2] = boxes[m + 2];
  boxes[k + 3] = boxes[m + 3];
  boxes[m] = a;
  boxes[m + 1] = b;
  boxes[m + 2] = c;
  boxes[m + 3] = d;

  const e = indices[i];
  indices[i] = indices[j];
  indices[j] = e;
}

// `hilbert`: compute hilbert codes.
//
// This is the function that takes a position in 2D space, `x` and `y`, and
// returns the hilbert value for that position.
//
// Umm yeah sorry I can't say anything else about this... it's black magic.
//
// Refer to the [C++ source](https://github.com/rawrunprotected/hilbert_curves)
// and the [original blog post](http://threadlocalmutex.com/?p=126) for any hope
// of understanding what's going on here!
/**
 * Fast Hilbert curve algorithm by http://threadlocalmutex.com/
 * Ported from C++ https://github.com/rawrunprotected/hilbert_curves (public domain)
 * @param {number} x
 * @param {number} y
 */
function hilbert(x, y) {
  let a = x ^ y;
  let b = 0xffff ^ a;
  let c = 0xffff ^ (x | y);
  let d = x & (y ^ 0xffff);

  let A = a | (b >> 1);
  let B = (a >> 1) ^ a;
  let C = (c >> 1) ^ (b & (d >> 1)) ^ c;
  let D = (a & (c >> 1)) ^ (d >> 1) ^ d;

  a = A;
  b = B;
  c = C;
  d = D;
  A = (a & (a >> 2)) ^ (b & (b >> 2));
  B = (a & (b >> 2)) ^ (b & ((a ^ b) >> 2));
  C ^= (a & (c >> 2)) ^ (b & (d >> 2));
  D ^= (b & (c >> 2)) ^ ((a ^ b) & (d >> 2));

  a = A;
  b = B;
  c = C;
  d = D;
  A = (a & (a >> 4)) ^ (b & (b >> 4));
  B = (a & (b >> 4)) ^ (b & ((a ^ b) >> 4));
  C ^= (a & (c >> 4)) ^ (b & (d >> 4));
  D ^= (b & (c >> 4)) ^ ((a ^ b) & (d >> 4));

  a = A;
  b = B;
  c = C;
  d = D;
  C ^= (a & (c >> 8)) ^ (b & (d >> 8));
  D ^= (b & (c >> 8)) ^ ((a ^ b) & (d >> 8));

  a = C ^ (C >> 1);
  b = D ^ (D >> 1);

  let i0 = x ^ y;
  let i1 = b | (0xffff ^ (i0 | a));

  i0 = (i0 | (i0 << 8)) & 0x00ff00ff;
  i0 = (i0 | (i0 << 4)) & 0x0f0f0f0f;
  i0 = (i0 | (i0 << 2)) & 0x33333333;
  i0 = (i0 | (i0 << 1)) & 0x55555555;

  i1 = (i1 | (i1 << 8)) & 0x00ff00ff;
  i1 = (i1 | (i1 << 4)) & 0x0f0f0f0f;
  i1 = (i1 | (i1 << 2)) & 0x33333333;
  i1 = (i1 | (i1 << 1)) & 0x55555555;

  return ((i1 << 1) | i0) >>> 0;
}
