# Magneto-Optical-Microsope

This repository contains all code, data structures, LaTeX summaries, notes and experimental logs for my internship at the FELIX Laboratory (Radboud University).  
The project focuses on magneto-optical microscopy and related experimental and theoretical work.

## Repository Structure:

```bash
Magneto-Optical-Microsope/
│
├── Code/                       # scripts, notebooks, simulatios
│   ├── Analysis/               # Experiment-specifieke analyses
│   ├── Simulations/            # Theoretical or Simulation code
│   │     ├── src/              # srs folder with all the files that contain the functions used in the simulation
│   │     │    ├── helper.py
│   │     │    ├── laser.py
│   │     │    ├── materials.py
│   │     │    ├── model_2tm.py
│   │     │    ├── model_llb.py
│   │     │    ├── postprocessing.py
│   │     │    └── simulation.py
│   │     │
│   │     └── main.py
│   │         
│   ├── README.md               # Explanation of code folder structure and dependencies
│   └── requirements.txt        # pip dependencies for venv
│
├── Experiments/                # data, notes, figures of experiments
│   ├── Figures/                # Schematics or important visuals
│   ├── Notes/                  # General notes, observations, settings, problems
│   ├── Processed_Data/         # Output of scripts in Code/
│   └── Raw_Data/               # Raw experimental data
│
├── Notes/                      # Flexible workspace / notebook
│   ├── Ideas.md                # Brainstorming / experimental ideas
│   ├── Meetings.md             # Notes from meetings
│   ├── Questions.md            # List of questions for supervisor
│   └── TODO.md                 # Todo list
│
├── Report/                     # Report and drafts
│   ├── Figures/                # Figures used in report
│   └── report.tex              # LateX file for final report
│
├── Summaries/                  # all summaries of papers
│   ├──bib/ 
│   │   └── bibliography.bib    # Bibliography with all citations
│   │
│   ├── PDFs/                   # Important papers(for easy acces)
│   └── papers_summaries.tex    # Summary of the papers
│
├── .gitignore                  # gitignore (venv, LateX)
└── README.md                   # General explanation about the repository and project
```

## Purpose of This Project

The goal of this internship is to develop experimental and analytical techniques for magneto-optical microscopy.  
This includes:

## Setup Virtual Environment

1. **Create virtual environment:**
   ```bash
   python -m venv venv
   ```
2. **Activate it:**
   ```bash
   source venv/bin/activate  # macOS/Linux
   venv\Scripts\activate     # Windows
   ```
3. **Install dependencies:**
   ```bash
   pip install -r Code/requirements.txt
   ```

## Repository Explanation

- **Summaries/** contains short LaTeX writeups of every paper I read  
- **Report/** will contain the internship final report (Latex)  
- **Code/** contains Python code, scripts, and utilities  
- **Experiments/** stores raw data, logs, setup parameters, processed data  
- **Notes/** includes planning, to-do lists, brainstorming, meeting notes  