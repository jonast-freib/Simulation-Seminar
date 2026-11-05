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

The test scenes are defined in `src/Scenes.cpp`. The app starts with the **dam break**: a fluid column in the lower-left corner of a tank collapses. The **hydrostatic column** fills the bottom of a narrower tank with fluid at rest; density should stay within about 1% of the rest density and pressure should grow linearly with depth.

Each frame advances 1/60 s of simulated time. The solver splits it into substeps whose size is chosen from the acoustic CFL, force and viscosity limits. Once per second the console prints the substep count, the current dt, the maximum density error and the Mach number.

## Runtime controls

The window displays debug sliders for:

- Rest density
- Speed of sound
- Viscosity
- Smoothing radius
- Gravity

Slider values are applied once per second. The buttons are:

- **Reset** rebuilds the active scene with its default parameters.
- **Dam** loads the dam-break scene.
- **Column** loads the hydrostatic-column scene.
- **Color** toggles the fluid colouring between the pressure map and plain blue (wall particles stay grey).

## Clean rebuild

To remove generated build output without touching source files, use a separate build directory or remove the local `build/` directory and configure again:

```sh
cmake -S . -B build-clean
cmake --build build-clean
```
