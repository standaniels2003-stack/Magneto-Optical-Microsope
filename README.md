# Magneto-Optical-Microsope

This repository contains all code, data structures, LaTeX summaries, notes and experimental logs for my internship at the FELIX Laboratory (Radboud University).  
The project focuses on magneto-optical microscopy and related experimental and theoretical work.

## Repository Structure:

```bash
Magneto-Optical-Microsope/
│
├── Code/ # scripts, notebooks, simulatios
│   ├── Analysis/ # Experiment-specifieke analyses
│   ├── Simulations/ # Theoretical or Simulation code
│   ├── src/ # Core helper modules / reusable functions
│   ├── README.md # Explanation of code folder structure and dependencies
│   └── requirements.txt # pip dependencies for venv
│
├── Experiments/ # data, notes, figures of experiments
│   ├── Figures/ # Schematics or important visuals
│   ├── Notes/ # General notes, observations, settings, problems
│   ├── Processed_Data/ # Output of scripts in Code/
│   └── Raw_Data/ # Raw experimental data
│
├── Notes/ # Flexible workspace / notebook
│   ├── Ideas.md # Brainstorming / experimental ideas
│   ├── Meetings.md # Notes from meetings
│   ├── Questions.md # List of questions for supervisor
│   └── TODO.md # Todo list
│
├── Report/ # Report and drafts
│   ├── Figures/ # Figures used in report
│   └── report.tex # LateX file for final report
│
├── Summaries/ # all summaries of papers
│   ├──bib/ 
│   │   └── bibliography.bib # Bibliography with all citations
│   ├── PDFs/ # Important papers(for easy acces)
│   └── papers_summaries.tex # Summary of the papers
│
├── .gitignore # gitignore (venv, LateX)
└── README.md # General explanation about the repository and project
```

## Purpose of This Project

The goal of this internship is to develop experimental and analytical techniques for magneto-optical microscopy.  
This includes:

## Setup Virtual Environment

1. Create virtual environment:
   python -m venv venv
2. Activate it:
   source venv/bin/activate  # macOS/Linux
   venv\Scripts\activate     # Windows
3. Install dependencies:
   pip install -r Code/requirements.txt