---
author: Lucas Timotheo Sanches
title: CapyrX
subtitle: Implementing the `Thornburg06` coordinate system
date: July 29, 2024
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

# Wave Toy