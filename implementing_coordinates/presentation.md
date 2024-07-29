---
author: Lucas Timotheo Sanches and Erik Schnetter
title: CapyrX
subtitle: Implementing the `Thornburg06` coordinate system
date: July 29, 2024

revealjs-url: "https://unpkg.com/reveal.js"
theme: dracula
controls: false
progress: true
slideNumber: true
width: 960
height: 700
margin: 0.0
---

# What's CapyrX ?

---

* `CapyrX` is a multipatch infrastructure for `CarpetX`.
* The name is a fusion of Capybara and Beaver (me and Erik).
* Today, we will be looking at how to implement a new coordinate system in it: The `Thornburg06` system

# The Repo:

---

The `CapyrX` code is an arrangement of thorns hosted in GitHub [here](https://github.com/lucass-carneiro/CapyrX)

---

To make the `Cactus` build system pull it and compile it automatically via the `GetComponents` + `simfactory` approach, simply add the following to your thorn list:

```bash
# CapyrX
!TARGET   = $ARR
!TYPE     = git
!URL      = https://github.com/lucass-carneiro/CapyrX
!REPO_PATH= $2
!CHECKOUT =
CapyrX/MultiPatch        # The meat
CapyrX/MultiPatchWaveToy # The test suite
CapyrX/TestMultiPatch    # The obligatory wave toy example
```

I will assume that you already did, if you plan on following along

# `Thornburg06` coordinates

## The Idea

* Cover the entire surface of the sphere without adding coordinate singularities.
* We do that by using 6 (hence the name) coordinate patches.
* Visualization: Imagine you have a malleable cube, and you blow air in it from the center in all directions. Each cube face becomes a spherical section. Together, these sections cover the whole $S^2$

## The Coordinates

1. We will be using `Mathematica`.
2. Go to `Capyrx/MultiPatch/src` and create a folder called `Thornburg06`
3. Lets create a subfolder called `resources` to keep our `Mathematica` files.
4. We will be looking at is `Thornburg06_deductions.nb` and `Thornburg06_visualize.nb`

# Matheamtica time!

# Include File

* Create `Capyrx/MultiPatch/src/Thornburg06/thornburg06.hxx`
* This is going to contain basic definitions and prototypes of functions we will need to fill later.

# Source and makefile

* Create `Capyrx/MultiPatch/src/Thornburg06/thornburg06.cxx`
* Edit `Capyrx/MultiPatch/src/make.code.defn`

---

```makefile
# Main make.code.defn file for thorn MultiPatch

# Source files in this directory
SRCS =						    \
	interpolate.cxx				\
	multipatch.cxx

# Subdirectories containing source files
SUBDIRS =        \
	cake         \
	cartesian    \
	cubed_sphere \
	swirl        \
	two_cubes    \
	Thornburg06  # <-- New!
```

---

* Create `Capyrx/MultiPatch/src/Thornburg06/make.code.defn`

```makefile
# make.code.defn file for the Thornburg06 patch system definition and tests

# Source files in this directory
SRCS = thornburg06.cxx tests.cxx
```

# Parameters

The `Thornburg06` coordinates will need the following parameters:

1. $r_0$: The inner boundary radius.
2. $r_1$: The outer boundary radius.
3. $n_r$: The number of radial cells.
4. $n_\theta$ and $n_\phi$: The number of angular cells

---

* We begin by adding these parameters in `MultiPatch/param.ccl`
* Next, we add the coordinate name `Thornburg06` as a possible parameter to the `patch_system` parameter declaration

---

* To access the `param.ccl` parameters, we need to edit the `PatchTransformations` in `MultiPatch/multipatch.hxx`
* We need to add a prototype for the `PatchSystem SetupThornburg06()` function.
* We need to edit the `PatchTransformations` constructor in `MultiPatch/multipatch.cxx`
* We need to add our coordinate system to `MultiPatch_Setup()`

# Coordinates source files

* To do this, we will need more Mathematica!
* Let us look at:
  1. `Thornburg06_transforms_codegen.nb`
  2. `Thornburg06_jacobians_codegen.nb`
* With these we can implement the patch system.

# Unit tests

* It is important to test the implementation of out coordinate system.
* In order to do that, we implement unit tests for each patch system.
* `CapyrX` has a set of common unit testing tools that can be found in the `MultiPatch/src/tests.hxx`

---

* Each coordinate further implements their unit tests in a `tests.cxx` folder. Users are responsible for implementing tests they think are relevant, but we have a few suggestions of useful tests that we have implemented in our patches

1. Given a random global point, test if the `get_owner_patch` function returns the correct patch
2. Given a random global point, test that `local2global(global2local(global)) == global`
3. Given a random local point and patch, test that `global2local(local2global(local, patch)) == (local, patch)`
4. Test that Jacobians and Jacobian derivatives obtained using finite differences are compatible with analytic implemented results.

---

* To schedule tests for running, we create a `C` linkage function called `MultiPatch_run_thornburg06_tests`
* We schedule this function in the `shcedule.ccl` file if the `run_tests` parameter is true. This is similar in all implemented patch systems.
* Since tests use random values, the parameter `test_repetitions` controls how many times to repeat each test, thus selecting new test values.
* To ensure reproducibility of the random numbers generated, the random seed and random generating algorithm are chosen at compile time.

---

* To actually run unit tests, we recommend the usage of the parameter files located in the `TestMultiPatch` thorn.
* While this is not strictly necessary, the `TestMultiPatch` thorn does additional testing of the system by performing interpolation and synchronization tests.
* In principle, tough, unit tests can be run independently of this thorn in production runs by simply setting the parameter to true.

# Wave Toy

* Let's run a wave toy using our new coordinate
* To do that, we have implemented the `MultiPatchWaveToy` thorn.
* This thorn does all the boring projections on the finite difference derivatives that are required when using `MultiPatch`.

---

* This serves both as test and example implementation for future users.
* We plan to eventually provide something akin to `GlobalDerivatives` to facilitate usage.
* The projections were calculated with the help of `Mathematica` once again.
* The code generation routines for derivative projections can be found in `MultiPatchWaveToy/resources/derivatives.nb`

# Plots

* Plotting `MultiPatch` simulation results can be painful.
* `Llama` had `VisIt` support and recently `Kuibit` support. For `CapyrX` we have no such luck.
* To ease the pain, `CapyrX` provides the `mpx` python script.
* I highly recommend using this script instead of rolling your own. I also recommend using this with `openPMD` data files (the `ASCII` support is not as developed).

---

The first step is installing the script dependencies. Navigate to the `CapyrX` repository folder and type

```bash
python -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

---

Let's dissect the common usage:

```bash
python3 mpx/mpx.py \
  plot-openpmd-slice \
  ~/simulations/mp/output-0000/sw_cake_opmd/sw_cake_opmd \
  multipatchwavetoy \
  state \
  u \
  0 \
  "[['01', True], ['02', True], ['03', True], ['04', True], ['00', True]]" \
  --verbose \
  --diverging \
  --varmin=-1.0 \
  --varmax=1.0 \
  --slice-coord=z \
  --slice-val=0.0 \
  --save
```
