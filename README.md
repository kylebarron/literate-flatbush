# literate-flatbush

At commit `4ab68d7e11607dba0e814efe50fb3583c69f90d6`:

Note: have a published "literate" version at multiple git tags, as well as a "main/index.html". At the top of the file, have a "versions" section that references/links to each of these commits.



# Overview

Before diving into the code, I want to provide a high level overview of how this R-Tree implementation works, so that hopefully the code will be easier to follow.

Flatbush is a really interesting R-Tree implementation because:

- Elegant, concise algorithm. Only 277 lines of code, excluding comments and blank lines, according to `cloc`.
- used as the basis for other projects like FlatGeobuf
- really fast
- really compact and memory efficient
- Only two buffer allocations: one for the main data buffer and a second intermediate one for the hilbert values during `finish`. In a lower-level language, this would mean _only two heap allocations_.
- it's a single contiguous buffer, which means that the index can be freely transfered from the main thread to a web worker without copies or serialization.
- The memory layout is defined exclusively by the number of items and the node size. So you can expect and plan for the amount of memory used.

There are only three parts to the data buffer:

- Header. This is just 8 bytes! One byte for the "magic" signifying this is a buffer containing a Flatbush index. The second byte is split, with four bits for the version and four bits for the `boxes`'s array type. Then two bytes for the node size and four for the number of items contained within.
- Boxes. This contains all the bounding box data, both for each _input geometry_ (tree leaf) and also for every _tree node_.
- Indices. Contains an ordering of `boxes`. This allows for retrieving the input ordering given some search query.

## Terminology

### Node size

The number of items contained within each node.

## Boxes

array data

## Indices

Note that the indices array is the length of `numNodes`, not `numItems`.

Therefore, the indices array serves two purposes. For positions 0-`numItems`, the indices array gives the original insertion order of the item that is located at that position.

For positions from `numItems` to `numNodes`, the indices array provides the structure of the tree.

## levelBounds

go from _leaves_ up

## bit shifts

There are a lot of bit shifts in the code, and that really confused me the first time I read it.

Bit shifts can be a useful performance optimization because they signify to the JavaScript JIT compiler that the output will _always be an integer_. If you do `Math.floor()`

// <img width="100%" src="hero.jpg" />
