# Testing Notes

The project currently uses a small CMake/CTest test executable for neighbor-search behavior.

## Run the suite

After configuring and building:

```sh
ctest --test-dir build --output-on-failure
```

Run the executable directly when debugging its diagnostic output:

```sh
./build/neighbor_search_test
```

## Covered behavior

The current tests exercise:

- Empty particle collections.
- A single particle with no neighbors.
- Particles inside the search radius.
- Particles outside the search radius.
- Changing the search radius and rebuilding the lists.
- A regular grid, including center, edge, and corner neighbor counts.

## Adding a test

Keep simulation-independent checks in `tests/` and register new executables or test commands in `CMakeLists.txt`. Tests should verify observable behavior rather than implementation details such as the internal storage layout of neighbor lists.

For geometry-based cases, state the radius and particle spacing near the fixture so the expected neighbor count is easy to audit. Rebuild the relevant target before running CTest when source or test files have changed.
