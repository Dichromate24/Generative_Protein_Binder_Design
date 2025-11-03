🧬 Generative Protein Binder Design Pipeline

This repository contains an end-to-end Python workflow for AI-assisted protein binder design using NVIDIA BioNeMo cloud models. The pipeline integrates multiple NVIDIA inference microservices (NIMs) — AlphaFold2, OpenFold3, ESMFold, RFdiffusion, ProteinMPNN, and AlphaFold-Multimer — to predict, design, and evaluate protein–protein interactions.


The workflow starts from a user-provided amino acid sequence (target protein):

1) Structure Prediction – User chosen model (AlphaFold2, OpenFold3, or ESMFold) predicts the 3D structure of the initial target protein.

2) Backbone Generation – RFdiffusion generates possible binder backbone structures conditioned on the predicted target structure.

3) Sequence Design – ProteinMPNN designs optimized amino acid sequences corresponding to the generated binder backbones.

4) Complex Evaluation – AlphaFold-Multimer predicts and assess potential binder–target complexes for structural compatibility and interaction quality.

5) Confidence Scoring – Calculate average pLDDT confidence scores for each predicted multimer model.
   

All steps are executed programmatically through REST API calls to NVIDIA NIM endpoints, with automatic polling, PDB file handling, py3Dmol visualization, and structured output saved to organized folders. The code is modular, human-readable, and easily adaptable for new protein design tasks or alternative models.
