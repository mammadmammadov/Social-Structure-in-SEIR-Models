# Social Structure in SEIR Models 🦠

![Built with R](https://img.shields.io/badge/Built%20with-R-276DC3?style=for-the-badge&logo=r)

This repository contains an **R** implementation of an SEIR (Susceptible, Exposed, Infected, Recovered) epidemic model.

Standard SEIR models often assume "random mixing," where individuals interact randomly like molecules of gas. This project refines that assumption by introducing a more realistic social structure, modeling how epidemics spread through households and regular contact networks.

## The Model: Core Concepts

This simulation models the spread of a disease by accounting for three distinct pathways of infection:

1.  **Households** 👨‍👩‍👧‍👦
    * Individuals are grouped into households of varying sizes (from 1 to `h_max`).
    * There is a high daily probability, $\alpha_{h}$, of infecting a household member.

2.  **Contact Networks** 🕸️
    * Beyond households, individuals have a regular, fixed network of contacts (e.g., friends, colleagues).
    * This network is randomly generated based on a "sociability" parameter ($\beta_i$) for each person.
    * Infection probability within this regular network is $\alpha_{c}$.

3.  **Random Mixing** 💨
    * A background level of random infection, $\alpha_{r}$, still occurs, representing general, non-regular interactions in the wider population.

## Key Functions 📊

The simulation is built around several core functions written in base R:

* **`h` vector generation**: A function that efficiently assigns all `n` individuals to households.
* **`get.net(beta, h, nc)`**: Generates the random contact network. It takes individual sociability parameters (`beta`) and household assignments (`h`) to create a unique contact list for each person, excluding their household members.
* **`nseir(...)`**: The main simulation engine. This function runs the SEIR model for `nt` days, processing new infections daily based on the three transmission pathways (household, network, and random).
* **Plotting Function**: A custom function to plot the S, E, I, and R population dynamics over time, as returned by `nseir`.

## Analysis: Investigating Social Structure

The primary script uses these functions to compare four distinct scenarios, illustrating the impact of social structure on epidemic dynamics.

The four scenarios compared are:

1.  **Full Model:** The complete model with all structures active (households, networks, and random mixing) and variable sociability (`beta`) for each person.
2.  **Random Mixing Only:** Simulates a traditional model by removing household and network effects ($\alpha_{h}=0$, $\alpha_{c}=0$) and setting $\alpha_{r}=0.04$.
3.  **Full Model (Constant $\beta$):** The full model, but with individual sociability removed (the `beta` parameter is set to a constant average for all individuals).
4.  **Random Mixing (Constant $\beta$):** The simplest model, combining both random-only mixing and constant `beta`.

The goal of this comparison is to visually and analytically determine how household and network structures, versus simple random mixing, affect the speed, peak, and overall size of an epidemic.

## 🚀 How to Run the Script (in VS Code)

### 🧩 Prerequisites

- **R:** You must have R installed on your system.  
- **Visual Studio Code:** The code editor.  
- **R Extension for VS Code:** Install the official R extension from the VS Code Marketplace.  

---

### 📥 Clone the Repository

First, clone the project repository to your local machine. This will provide you with both the R script (`proj1.r`) and the required `shakespeare.txt` file.

```bash
https://github.com/mammadmammadov/Social-Structure-in-SEIR-Models.git
```

### 💻 Open in VS Code

Navigate to the project folder and open it in Visual Studio Code:

```bash
cd Social-Structure-in-SEIR-Models
code .
```

### ▶️ Run the Script

1. Open the `proj2.r` file in the editor.
2. With the file open, you can run the entire script by pressing `Ctrl + Shift + S`

