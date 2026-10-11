# Building and Running

## Configure

Configure an out-of-source build from the repository root:

```sh
cmake -S . -B build
```

CMake locates SFML through its package configuration. If SFML is installed in a non-standard prefix, provide that prefix during configuration:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/sfml
```

## Compile

Build both the application and test executable:

```sh
cmake --build build
```

For an incremental build of the simulator only:

```sh
cmake --build build --target sph_solver
```

## Launch

Run the simulator from the repository root so relative runtime resources, if added later, resolve predictably:

```sh
./build/sph_solver
```

The current scene is created in `src/main.cpp`. The active startup preset is a single movable particle above a fixed U-shaped block.

## Runtime controls

The window displays sliders for:

- Rest density
- Gas constant
- Viscosity
- Smoothing radius
- Timestep
- Gravity

The **Reset** button restores particle positions and clears velocity and acceleration for the active scene. Fixed particles remain fixed by the solver.

## Clean rebuild

To remove generated build output without touching source files, use a separate build directory or remove the local `build/` directory and configure again:

```sh
cmake -S . -B build-clean
cmake --build build-clean
```
