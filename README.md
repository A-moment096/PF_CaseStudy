# PF_CaseStudy: Phase Field Simulation Example Code

Case studies from *Programming Phase-Field Modeling*, but with C++.

## Summary

This repository presents the C++ code for the first 5 case studies in book *Programming Phase-Field Modeling*, by S. Bulent Biner.

This implementation is based on object orientated programming paradigm, with slightly little C++ templates.

Please take the code presented here as an exercise and attempt for exploring the possibility for a generalized phase field program.

## Code Overview

There are six folders, with five folders containing case studies from 1 to 5, and a folder named **include** containing all headers used in each case study.

### Case Studies

Below are the details of each case study. Please notice that several case study is not fully implemented as the example in the book.

- CaseStudy_1: A-B alloy spinodal decomposition;
    - Cahn-Hilliard (C-H) equation example;
    - Double well bulk energy;
    - Simple gradient interface energy;
- CaseStudy_2: Grain growth (Partially implemented with only two grains growth);
    - Allen-Cahn (A-C) equation example;
    - Polynomial bulk energy for multiple grain/phase variables: half-double-well energy for single phase variable added with interaction from other phase variables;
    - Simple gradient interface energy;
    
- CaseStudy_3: Solid state sintering
    - C-H equation for density and A-C equation for phase variable;
    - Double well energy from density;
    - Polynomial energy from phase variable coupled with density;
    - Simple gradient interface energy for both density and phase variable;
    - Complex diffusion, using interpolation function;
- CaseStudy_4: Dendrite growth
    - A-C equation for phase variable;
    - Conservation law of enthalpy for temperature evolution, with simple latent heat contribution;
    - Double well bulk energy for phase variable coupled with supercooling;
    - Gradient interface energy with anisotrophic gradient coefficient;
- CaseStudy_5: Cell motility
    - Soft cell motion;
    - A-C equation coupled with velocity field;
    - Velocity field decomposed into surrounding action part and active part;
    - Complex bulk energy coupled with constant volume restriction;
    - Interface energy describing interaction between cells (part of bulk energy form CaseStudy_2);

### Class framework

The code of these five case studies are all based on a set of class and method. Here is a breif introduction to the framework.

**Overall Structure** 

The whole data structure is based on Mesh - MeshNode - EntryNode - Entry model. It means that:

- There should be only one mesh with several points;
- On each point there will be a mesh node;
- In each mesh node, there should be several entry nodes representing different physical properties, such as phase variable, density, concentration.
- In each entry node, there are several entries for different variables.

For example, if there is a simulation of two concentration and one phase variable, with simulation region size of $100\times 100$, then there should be exactly one mesh, 10000 mesh node; in each mesh node there are two active entry node, for concentration and phase variable, and for concentration node, there will be two entry, element 1 and element 2; as for phase variable, there will be only one entry.

**File Structure**

With such structure design, the file structure is quite direct and simple.

- BaseEntry.hh: All kinds of entry will be here, derived from a parent class `BaseEntry` defined here.
- BaseNode.hh: All kinds of node defined, also derived from a parent class `BaseNode` defined here as well.
- MeshNode.hh: A mesh node to manage all base node. Can be considered as a *node manager*
- PFM.hh: A header for including everything that will be used in each case study.
- PFMTools.hh: A header for several utilities, such as timer or boundary value check, which are less relative to the data structure described before.
- SimulationMesh.hh: A header to define the whole mesh.

And, with such data structure, algorithms related to each structure become a part of method for each class. For instance, the algorithm to calculate laplacian is method of `SimulationMesh`, while update value of each node will be traced back to the method of `BaseNode` since all properties should be updated.
