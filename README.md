🧬 Generative Protein Binder Design Pipeline

This repository contains an end-to-end Python workflow for AI-assisted protein binder design using NVIDIA BioNeMo cloud models. The pipeline integrates multiple neural inference microservices (NIMs) — AlphaFold2, OpenFold3, ESMFold, RFdiffusion, ProteinMPNN, and AlphaFold-Multimer — to predict, design, and evaluate protein–protein interactions.

The workflow starts from a user-provided amino acid sequence (target protein):

1) Predicts the 3D structure of the input target protein using a selected model (AlphaFold2, OpenFold3, ESMFold).

2) Generates new binder candidates with RFdiffusion based on the predicted structure.

3) Designs optimized sequences for the binders via ProteinMPNN.

4) Evaluates potential complexes through AlphaFold-Multimer for structural compatibility.

5) Computes average pLDDT confidence scores for each predicted multimer.



All steps are executed programmatically through REST API calls to NVIDIA NIM endpoints, with automatic polling, PDB file handling, py3Dmol visualization, and structured output saved to organized folders. The code is modular, human-readable, and easily adaptable for new protein design tasks or alternative models.
