# AlphaFold ANGPTL-Lipase Interaction Pipeline

**Author:** Alexander Dumas
**Last Updated:** July 2025

## Overview
This repository contains scripts to analyze AlphaFold3 predictions of Angiopoietin like proteins and their ability to interact with each other and members of the pancreatic lipase family

Key features include:
- Automated AlphaFold submission using Playwright (semi-automated for heterotrimer analysis)
- Extraction of interacting residues using PyMOL and distance-based filtering.
- Bar plot generation of residue-level interaction frequencies across multiple models.
- Highlighting of conserved or catalytic residues across multiple lipases.
- Jaccard similarity calculations to assess prediction consistency.

## Repository structure
