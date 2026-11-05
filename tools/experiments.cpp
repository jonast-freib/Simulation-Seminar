// Headless experiment runner for the seminar paper.
//
// Runs the dam-break and hydrostatic-column scenes for the default setup and
// for parameter sweeps (speed of sound, viscosity, h/dx) and writes CSV files
// that Latex/plots/make_plots.py turns into figures.
//
// Usage: sph_experiments <output-directory>
//
// Coordinates in all snapshot files are relative to the tank: x is measured
// from the left wall surface, y is the height above the floor.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "Scenes.h"
#include "SPHSolver.h"

namespace
{
constexpr double windowWidth = 1700.0;
constexpr double windowHeight = 900.0;
constexpr double frameTime = 1.0 / 60.0;

std::filesystem::path outputDirectory;
std::mutex logMutex;

void log(const std::string &message)
{
    std::lock_guard<std::mutex> lock(logMutex);
    std::printf("%s\n", message.c_str());
    std::fflush(stdout);
}

std::ofstream openCsv(const std::string &name, const std::string &header)
{
    std::ofstream file(outputDirectory / name);
    file << header << '\n';
    return file;
}

std::string formatValue(double value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

// A scene together with a solver configured for it.
struct Simulation
{
    Scene scene;
    SPHSolver solver;
    double time = 0.0;
    long substeps = 0;
    long frames = 0;

    explicit Simulation(Scene s) : scene(std::move(s))
    {
        applySceneParameters(solver, scene.parameters);
        solver.setBounds(windowWidth, windowHeight);
        // One zero-length step so density and pressure are defined at t = 0.
        solver.step(scene.particles, 0.0);
    }

    void advanceFrame()
    {
        substeps += solver.advance(scene.particles, frameTime);
        ++frames;
        time += frameTime;
    }

    double floorY() const { return scene.tank.top + scene.tank.height; }
    double heightAboveFloor(const Particle &p) const { return floorY() - p.position.y; }
    double xFromLeftWall(const Particle &p) const { return p.position.x - scene.tank.left; }
    double rho0() const { return scene.parameters.restDensity; }
    double gravity() const { return scene.parameters.gravity; }
    // Hydrostatic pressure scale at the bottom of the initial fluid body.
    double pressureScale() const { return rho0() * gravity() * scene.fluidHeight; }
};

struct DensityStats
{
    double maxError = 0.0;  // max (rho / rho0 - 1) over fluid particles
    double meanError = 0.0; // mean (rho / rho0 - 1)
    double meanKinetic = 0.0;
    double maxSpeed = 0.0;
};

DensityStats densityStats(const Simulation &sim)
{
    DensityStats stats;
    int count = 0;
    for (const Particle &p : sim.scene.particles)
    {
        if (p.fixed)
            continue;
        const double error = p.density / sim.rho0() - 1.0;
        const double speed2 = p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y;
        stats.maxError = std::max(stats.maxError, error);
        stats.meanError += error;
        stats.meanKinetic += 0.5 * speed2;
        stats.maxSpeed = std::max(stats.maxSpeed, std::sqrt(speed2));
        ++count;
    }
    if (count > 0)
    {
        stats.meanError /= count;
        stats.meanKinetic /= count;
    }
    return stats;
}

// Leading edge of the collapsing column along the floor, measured from the
// left wall. Only particles within 3 spacings of the floor count, so splashes
// do not define the front. Half a spacing is added for the particle extent,
// so the initial front equals the column width a.
double damFront(const Simulation &sim)
{
    const double dx = sim.solver.getParticleSpacing();
    double front = 0.0;
    for (const Particle &p : sim.scene.particles)
    {
        if (p.fixed || sim.heightAboveFloor(p) > 3.0 * dx)
            continue;
        front = std::max(front, sim.xFromLeftWall(p) + 0.5 * dx);
    }
    return front;
}

void writeSnapshot(const Simulation &sim, const std::string &name,
                   const std::vector<double> *smoothedPressure = nullptr)
{
    std::ofstream file = openCsv(name, "x,y,fixed,density_error,pressure,pressure_smoothed,speed");
    for (std::size_t i = 0; i < sim.scene.particles.size(); ++i)
    {
        const Particle &p = sim.scene.particles[i];
        const double smoothed = smoothedPressure ? (*smoothedPressure)[i] : p.pressure;
        file << sim.xFromLeftWall(p) << ',' << sim.heightAboveFloor(p) << ',' << (p.fixed ? 1 : 0) << ','
             << (p.fixed ? 0.0 : p.density / sim.rho0() - 1.0) << ',' << p.pressure << ',' << smoothed << ','
             << std::hypot(p.velocity.x, p.velocity.y) << '\n';
    }
}

// Exponential moving average of each particle's pressure (time constant tau),
// the same smoothing the renderer uses to suppress acoustic flicker.
void updateSmoothedPressure(const Simulation &sim, std::vector<double> &smoothed, double tau)
{
    const double blend = 1.0 - std::exp(-frameTime / tau);
    if (smoothed.size() != sim.scene.particles.size())
    {
        smoothed.clear();
        for (const Particle &p : sim.scene.particles)
            smoothed.push_back(p.pressure);
        return;
    }
    for (std::size_t i = 0; i < smoothed.size(); ++i)
        smoothed[i] += blend * (sim.scene.particles[i].pressure - smoothed[i]);
}

void writeSceneMeta()
{
    std::ofstream file = openCsv("scene_meta.csv", "scene,key,value");
    for (int which = 0; which < 2; ++which)
    {
        const Scene scene = which == 0 ? makeDamBreakScene(windowWidth, windowHeight)
                                       : makeHydrostaticScene(windowWidth, windowHeight);
        const std::string name = which == 0 ? "dam" : "hydro";
        int fluid = 0;
        int wall = 0;
        for (const Particle &p : scene.particles)
            (p.fixed ? wall : fluid)++;
        const SceneParameters &q = scene.parameters;
        const std::vector<std::pair<std::string, double>> rows = {
            {"tank_width", scene.tank.width},
            {"tank_height", scene.tank.height},
            {"closed_top", scene.tank.closedTop ? 1.0 : 0.0},
            {"fluid_width", scene.fluidWidth},
            {"fluid_height", scene.fluidHeight},
            {"particle_spacing", q.smoothingRadius / q.kernelSpacingRatio},
            {"smoothing_length", q.smoothingRadius},
            {"kernel_spacing_ratio", q.kernelSpacingRatio},
            {"rest_density", q.restDensity},
            {"speed_of_sound", q.speedOfSound},
            {"viscosity", q.viscosity},
            {"gravity", q.gravity},
            {"fluid_particles", static_cast<double>(fluid)},
            {"wall_particles", static_cast<double>(wall)},
        };
        for (const auto &[key, value] : rows)
            file << name << ',' << key << ',' << value << '\n';
    }
}

// E1: hydrostatic column with the default setup.
void runHydrostatic()
{
    Simulation sim(makeHydrostaticScene(windowWidth, windowHeight));
    std::ofstream series = openCsv("hydro_timeseries.csv",
                                   "t,max_density_error,mean_density_error,mean_kinetic_energy,max_speed,substeps,dt");

    const double settleTime = 8.0;
    const double averageTime = 2.0;
    std::vector<double> pressureSum(sim.scene.particles.size(), 0.0);
    std::vector<double> heightSum(sim.scene.particles.size(), 0.0);
    int averagedFrames = 0;

    while (sim.time < settleTime + averageTime - 1e-9)
    {
        const long before = sim.substeps;
        sim.advanceFrame();
        const DensityStats stats = densityStats(sim);
        series << sim.time << ',' << stats.maxError << ',' << stats.meanError << ',' << stats.meanKinetic << ','
               << stats.maxSpeed << ',' << (sim.substeps - before) << ',' << sim.solver.getLastTimeStep() << '\n';

        if (sim.time > settleTime)
        {
            for (std::size_t i = 0; i < sim.scene.particles.size(); ++i)
            {
                pressureSum[i] += sim.scene.particles[i].pressure;
                heightSum[i] += sim.heightAboveFloor(sim.scene.particles[i]);
            }
            ++averagedFrames;
        }
    }

    // Depth below the initial free surface; the column compresses by < 1 %.
    std::ofstream profile = openCsv("hydro_profile.csv", "x,depth,pressure,pressure_time_averaged,hydrostatic");
    for (std::size_t i = 0; i < sim.scene.particles.size(); ++i)
    {
        const Particle &p = sim.scene.particles[i];
        if (p.fixed)
            continue;
        const double depth = sim.scene.fluidHeight - heightSum[i] / averagedFrames;
        profile << sim.xFromLeftWall(p) << ',' << depth << ',' << p.pressure << ','
                << pressureSum[i] / averagedFrames << ',' << sim.rho0() * sim.gravity() * depth << '\n';
    }
    writeSnapshot(sim, "hydro_snapshot.csv");
    log("E1 hydrostatic column done");
}

// E2: dam break with the default setup.
void runDamBreak()
{
    Simulation sim(makeDamBreakScene(windowWidth, windowHeight));
    const double a = sim.scene.fluidWidth;
    const double timeScale = std::sqrt(sim.gravity() / a); // T = t * sqrt(g / a)
    std::ofstream front = openCsv("dam_front.csv", "t,T,Z_over_a,max_density_error,mach");

    const std::vector<double> snapshotTimes = {0.0, 1.0, 2.0, 3.0, 5.0};
    std::size_t nextSnapshot = 0;
    std::vector<double> smoothed;
    const double endTime = 6.0 / timeScale;

    auto record = [&]()
    {
        const DensityStats stats = densityStats(sim);
        const double T = sim.time * timeScale;
        front << sim.time << ',' << T << ',' << damFront(sim) / a << ',' << stats.maxError << ','
              << stats.maxSpeed / sim.scene.parameters.speedOfSound << '\n';
        if (nextSnapshot < snapshotTimes.size() && T >= snapshotTimes[nextSnapshot] - 1e-9)
        {
            char name[64];
            std::snprintf(name, sizeof(name), "dam_snapshot_T%.0f.csv", snapshotTimes[nextSnapshot]);
            writeSnapshot(sim, name, &smoothed);
            ++nextSnapshot;
        }
    };

    updateSmoothedPressure(sim, smoothed, 0.05);
    record();
    while (sim.time < endTime - 1e-9)
    {
        sim.advanceFrame();
        updateSmoothedPressure(sim, smoothed, 0.05);
        record();
    }
    log("E2 dam break done");
}

// Results of one sweep configuration.
struct SweepResult
{
    std::string parameter;
    double value = 0.0;
    // dam break
    double damMaxDensityError = 0.0;  // max over the run
    double damMeanMaxDensityError = 0.0; // time average of the per-frame max
    double damSubstepsPerFrame = 0.0;
    double damRuntime = 0.0;          // wall-clock seconds per simulated second
    double damFrontAtT15 = 0.0;
    // hydrostatic column (after settling)
    double hydroMaxDensityError = 0.0;
    double hydroPressureDeviation = 0.0;        // mean |p - rho0 g d| / (rho0 g H), instantaneous
    double hydroAveragedPressureDeviation = 0.0; // same for 1 s time-averaged pressure
    double hydroFrameToFrameChange = 0.0;      // mean |p(t) - p(t - frame)| / (rho0 g H)
    double hydroMeanNeighbors = 0.0;
};

std::string sweepTag(const std::string &parameter, double value)
{
    return parameter + "_" + formatValue(value);
}

// E3/E4: one point of a parameter sweep. Runs a short dam break and a
// hydrostatic column and writes the density snapshots used for the
// small-multiples figure and the layer histogram for the pairing study.
SweepResult runSweepPoint(const std::string &parameter, double value, const SceneOptions &options)
{
    SweepResult result;
    result.parameter = parameter;
    result.value = value;
    const std::string tag = sweepTag(parameter, value);

    // --- dam break up to T = 3 ---
    {
        Simulation sim(makeDamBreakScene(windowWidth, windowHeight, options));
        const double a = sim.scene.fluidWidth;
        const double timeScale = std::sqrt(sim.gravity() / a);
        std::ofstream front = openCsv("sweep_front_" + tag + ".csv", "T,Z_over_a,max_density_error");
        bool snapshotWritten = false;
        bool frontRecorded = false;
        double errorSum = 0.0;
        const auto start = std::chrono::steady_clock::now();
        while (sim.time * timeScale < 3.0 - 1e-9)
        {
            sim.advanceFrame();
            const double T = sim.time * timeScale;
            const DensityStats stats = densityStats(sim);
            result.damMaxDensityError = std::max(result.damMaxDensityError, stats.maxError);
            errorSum += stats.maxError;
            front << T << ',' << damFront(sim) / a << ',' << stats.maxError << '\n';
            if (!frontRecorded && T >= 1.5 - 1e-9)
            {
                result.damFrontAtT15 = damFront(sim) / a;
                frontRecorded = true;
            }
            if (!snapshotWritten && T >= 2.0 - 1e-9)
            {
                writeSnapshot(sim, "density_dam_" + tag + ".csv");
                snapshotWritten = true;
            }
        }
        const double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        result.damMeanMaxDensityError = errorSum / sim.frames;
        result.damSubstepsPerFrame = static_cast<double>(sim.substeps) / sim.frames;
        result.damRuntime = elapsed / sim.time;
    }

    // --- hydrostatic column: settle 5 s, then measure for 2 s ---
    {
        Simulation sim(makeHydrostaticScene(windowWidth, windowHeight, options));
        const double settleTime = 5.0;
        const double measureTime = 2.0;
        while (sim.time < settleTime - 1e-9)
            sim.advanceFrame();

        const std::size_t n = sim.scene.particles.size();
        const double scale = sim.pressureScale();
        const double surfaceHeight = sim.scene.fluidHeight;
        std::vector<double> previous(n), pressureSum(n, 0.0), depthSum(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            previous[i] = sim.scene.particles[i].pressure;

        // Height histogram above the floor (1 px bins) in the central region,
        // away from the side walls: reveals layering / pairing.
        const int histogramBins = 80;
        std::vector<double> histogram(histogramBins, 0.0);
        const double margin = 4.0 * sim.solver.getParticleSpacing();

        double deviationSum = 0.0;
        double changeSum = 0.0;
        double neighborSum = 0.0;
        long samples = 0;
        int frames = 0;
        while (sim.time < settleTime + measureTime - 1e-9)
        {
            sim.advanceFrame();
            ++frames;
            for (std::size_t i = 0; i < n; ++i)
            {
                const Particle &p = sim.scene.particles[i];
                if (p.fixed)
                    continue;
                const double depth = surfaceHeight - sim.heightAboveFloor(p);
                deviationSum += std::abs(p.pressure - sim.rho0() * sim.gravity() * depth) / scale;
                changeSum += std::abs(p.pressure - previous[i]) / scale;
                neighborSum += static_cast<double>(sim.solver.getNeighborCount(i));
                pressureSum[i] += p.pressure;
                depthSum[i] += depth;
                previous[i] = p.pressure;
                ++samples;

                const double x = sim.xFromLeftWall(p);
                if (x > margin && x < sim.scene.tank.width - margin)
                {
                    const int bin = static_cast<int>(sim.heightAboveFloor(p));
                    if (bin >= 0 && bin < histogramBins)
                        histogram[bin] += 1.0;
                }
            }
            result.hydroMaxDensityError = std::max(result.hydroMaxDensityError, densityStats(sim).maxError);
        }

        double averagedDeviation = 0.0;
        int fluidCount = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (sim.scene.particles[i].fixed)
                continue;
            const double depth = depthSum[i] / frames;
            averagedDeviation += std::abs(pressureSum[i] / frames - sim.rho0() * sim.gravity() * depth) / scale;
            ++fluidCount;
        }
        result.hydroPressureDeviation = deviationSum / samples;
        result.hydroFrameToFrameChange = changeSum / samples;
        result.hydroAveragedPressureDeviation = averagedDeviation / fluidCount;
        result.hydroMeanNeighbors = neighborSum / samples;

        std::ofstream layers = openCsv("layers_" + tag + ".csv", "height,count_fraction");
        double total = 0.0;
        for (double c : histogram)
            total += c;
        for (int b = 0; b < histogramBins; ++b)
            layers << b + 0.5 << ',' << (total > 0.0 ? histogram[b] / total : 0.0) << '\n';
        writeSnapshot(sim, "density_hydro_" + tag + ".csv");
    }

    log("sweep " + tag + " done");
    return result;
}

void writeSweep(const std::string &name, const std::vector<SweepResult> &results)
{
    std::ofstream file = openCsv(
        name,
        "parameter,value,dam_max_density_error,dam_mean_max_density_error,dam_substeps_per_frame,"
        "dam_runtime_per_sim_second,dam_front_T1.5,hydro_max_density_error,hydro_pressure_deviation,"
        "hydro_averaged_pressure_deviation,hydro_frame_to_frame_change,hydro_mean_neighbors");
    for (const SweepResult &r : results)
    {
        file << r.parameter << ',' << r.value << ',' << r.damMaxDensityError << ',' << r.damMeanMaxDensityError
             << ',' << r.damSubstepsPerFrame << ',' << r.damRuntime << ',' << r.damFrontAtT15 << ','
             << r.hydroMaxDensityError << ',' << r.hydroPressureDeviation << ','
             << r.hydroAveragedPressureDeviation << ',' << r.hydroFrameToFrameChange << ','
             << r.hydroMeanNeighbors << '\n';
    }
}

// Run jobs on all hardware threads.
void runParallel(std::vector<std::function<void()>> &jobs)
{
    std::atomic<std::size_t> next{0};
    const unsigned threadCount = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    for (unsigned t = 0; t < threadCount; ++t)
    {
        threads.emplace_back([&]()
                             {
            for (std::size_t job = next++; job < jobs.size(); job = next++)
                jobs[job](); });
    }
    for (std::thread &thread : threads)
        thread.join();
}
} // namespace

int main(int argc, char **argv)
{
    outputDirectory = argc > 1 ? argv[1] : "Latex/data";
    std::filesystem::create_directories(outputDirectory);
    std::printf("Writing experiment data to %s\n", outputDirectory.c_str());

    const std::vector<double> speedOfSoundFactors = {2.5, 5.0, 10.0, 20.0};
    const std::vector<double> viscosities = {50.0, 150.0, 400.0, 1000.0};
    const std::vector<double> kernelRatios = {1.2, 1.3, 1.6, 2.0};

    std::vector<SweepResult> cResults(speedOfSoundFactors.size());
    std::vector<SweepResult> viscResults(viscosities.size());
    std::vector<SweepResult> ratioResults(kernelRatios.size());

    std::vector<std::function<void()>> jobs;
    jobs.push_back(writeSceneMeta);
    jobs.push_back(runHydrostatic);
    jobs.push_back(runDamBreak);
    // Most expensive jobs first so the thread pool finishes evenly.
    for (std::size_t i = 0; i < kernelRatios.size(); ++i)
    {
        jobs.push_back([&, i]()
                       {
            SceneOptions options;
            options.kernelSpacingRatio = kernelRatios[i];
            ratioResults[i] = runSweepPoint("hdx", kernelRatios[i], options); });
    }
    for (std::size_t i = 0; i < speedOfSoundFactors.size(); ++i)
    {
        jobs.push_back([&, i]()
                       {
            SceneOptions options;
            options.speedOfSoundFactor = speedOfSoundFactors[i];
            cResults[i] = runSweepPoint("c", speedOfSoundFactors[i], options); });
    }
    for (std::size_t i = 0; i < viscosities.size(); ++i)
    {
        jobs.push_back([&, i]()
                       {
            SceneOptions options;
            options.viscosity = viscosities[i];
            viscResults[i] = runSweepPoint("visc", viscosities[i], options); });
    }

    const auto start = std::chrono::steady_clock::now();
    runParallel(jobs);

    writeSweep("sweep_c.csv", cResults);
    writeSweep("sweep_visc.csv", viscResults);
    writeSweep("sweep_hdx.csv", ratioResults);

    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("All experiments finished in %.0f s\n", elapsed);
    return 0;
}
