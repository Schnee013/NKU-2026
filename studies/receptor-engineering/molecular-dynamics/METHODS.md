# Simulation and analysis methods

## Molecular systems

The study contains full-length Gpr2 TWF154 WT, N84V and V368R receptors, each with ascr18 or ascr3. Every receptor contains residues 1–406 and one ligand. Prepare full-length receptor predictions using the [AF3 workflow guide](../../../REPRODUCIBILITY.md#alphafold-3-predictions). The original six systems use their matched full-length receptor predictions, docked ligand poses and preparation; the two directed starts derive from the V368R–ascr3 parent described below. The comparisons describe these six starting conformations and two additional V368R–ascr3 starting poses. Ascarosides were modelled with a deprotonated carboxylate. The structural model includes a Cys99–Cys173 disulfide bond.

WT–ascr18 retains the native lid with a Trp174 side-chain rotation of χ1 +160° and χ2 −30°. For WT–ascr3 and V368R–ascr3, residues 164–178 were rebuilt while residues 163/179 and the remaining receptor heavy atoms were fixed during rebuilding. Each rebuilt system used its own receptor and ligand pose. The N84V arms and V368R–ascr18 retain their predicted lids.

The protein force field is CHARMM36m in the July 2022 CHARMM36 parameter distribution; ligand parameters use CGenFF 4.6, with CHARMM-GUI 1.7 topology conversion. The membrane contains POPC, the solvent is CHARMM TIP3P water, and counterions are K⁺ and Cl⁻. The molecule counts below were read from the production topologies.

| System | Total atoms | POPC | Water molecules | K⁺ | Cl⁻ |
|---|---:|---:|---:|---:|---:|
| WT–ascr18 | 145,776 | 180 | 38,327 | 104 | 117 |
| N84V–ascr18 | 146,032 | 177 | 38,545 | 105 | 118 |
| N84V–ascr3 | 146,021 | 177 | 38,544 | 105 | 118 |
| V368R–ascr18 | 145,979 | 177 | 38,525 | 105 | 119 |
| V368R–ascr3, rebuilt lid | 145,905 | 177 | 38,503 | 105 | 119 |
| WT–ascr3, rebuilt lid | 145,681 | 180 | 38,298 | 104 | 117 |
| V368R–ascr3, Arg368-directed | 145,881 | 177 | 38,495 | 105 | 119 |
| V368R–ascr3, Arg189-directed | 145,881 | 177 | 38,495 | 105 | 119 |

[System configurations](configs/simulation/systems.json) associate these compositions with the dataset IDs. [Computational environments](ENVIRONMENT.md) give the production and analysis hardware and software versions.

## Preparation and production

System preparation used steepest-descent energy minimization, heating, a constant-volume temperature hold and five constant-pressure stages. The [per-system configuration](configs/simulation/systems.json) records molecule counts, minimization file order and precision, and the heating, hold, NPT and production files.

Regular minimization targets a maximum force of 1,000 kJ mol⁻¹ nm⁻¹ within at most 100,000 steps. [em.mdp](configs/simulation/em.mdp) uses a maximum step of 0.001 nm; [em-small-step.mdp](configs/simulation/em-small-step.mdp) uses 0.0001 nm. WT–ascr18 and both rebuilt-lid systems use the former; N84V–ascr18, N84V–ascr3 and V368R–ascr18 use the latter.

N84V–ascr3 first undergoes [constrained minimization](configs/simulation/em-ligand-constrained.mdp): at most 10,000 steps, a 0.0001 nm step, and LINCS order 8 with four iterations. Its [ligand topology](configs/simulation/ascr3-minimization-constrained.itp) constrains ligand heavy-atom bonds. V368R–ascr18 first undergoes [frozen-ligand minimization](configs/simulation/em-ligand-frozen.mdp), also with at most 10,000 steps and a 0.0001 nm step. All 54 ligand atoms are frozen in three dimensions; [em-frozen-ligand.ndx](configs/simulation/em-frozen-ligand.ndx) identifies them with one-based atom indices. Both preliminary stages are followed by regular minimization with ordinary ligand topology and mobile ligand coordinates. Heating and subsequent stages use the ordinary ligand topology.

| Stage | Parameter file | Duration / ps | Step / fs | Preparation clock / ps |
|---|---|---:|---:|---|
| Heating | [heat.mdp](configs/simulation/heat.mdp) | 50 | 1 | 0–50 |
| Temperature hold | [hold.mdp](configs/simulation/hold.mdp) | 75 | 1 | 50–125 |
| NPT 1 | [npt01.mdp](configs/simulation/npt01.mdp) | 125 | 1 | 125–250 |
| NPT 2 | [npt02.mdp](configs/simulation/npt02.mdp) | 500 | 2 | 250–750 |
| NPT 3 | [npt03.mdp](configs/simulation/npt03.mdp) | 500 | 2 | 750–1,250 |
| NPT 4 | [npt04.mdp](configs/simulation/npt04.mdp) | 500 | 2 | 1,250–1,750 |
| NPT 5 | [npt05.mdp](configs/simulation/npt05.mdp) | 1,000 | 2 | 1,750–2,750 |

Heating increases all three temperature groups linearly from 50 to 303.15 K. The original six systems generate velocities at 50 K with `gen_seed = 8428368`; the Arg368-directed and Arg189-directed systems use fresh seeds `202610153` and `202610159`, respectively. Later stages continue their own velocities. NPT stages use 303.15 K and semiisotropic 1 bar C-rescale coupling.

Position-restraint coefficients are in kJ mol⁻¹ nm⁻². `DIHRES_FC` is in kJ mol⁻¹ rad⁻² for the topology's [restrained-dihedral terms](https://manual.gromacs.org/2025.4/reference-manual/functions/restraints.html#dihedral-restraints). The MDP `define` entries set the following coefficients; restrained pressure-coupled stages use `refcoord_scaling = all`.

| Stage | Receptor backbone | Receptor side chains | Lipid | Ligand | `DIHRES_FC` |
|---|---:|---:|---:|---:|---:|
| Minimization / heating | 4,000 | 2,000 | 1,000 | 1,000 | 1,000 |
| Hold | 1,000 | 500 | 400 | 100 | 400 |
| NPT 1 | 1,000 | 500 | 400 | 100 | 200 |
| NPT 2 | 500 | 200 | 200 | 50 | 200 |
| NPT 3 | 200 | 50 | 40 | 10 | 100 |
| NPT 4 | 50 | 0 | 0 | 0 | 0 |
| NPT 5 / production | 0 | 0 | 0 | 0 | 0 |

Selected energy minimization used double precision for the four native/predicted-lid systems and mixed precision for the two rebuilt-lid systems. The two directed-start systems also use mixed-precision minimization. All eight production simulations used mixed precision. Production continued from each system's equilibrated coordinates and velocities, with its production time reset to zero. Each system contributed one 100 ns trajectory, as recorded in [inputs.json](inputs.json).

| Parameter | Production setting |
|---|---|
| Integrator and time step | Leap-frog MD; 0.002 ps (2 fs) |
| Steps and duration | 50,000,000 steps; 100,000 ps |
| Periodic boundaries | Three-dimensional, `pbc = xyz` |
| Temperature coupling | Stochastic velocity rescaling, 303.15 K; independent SOLU (receptor + ligand), MEMB (POPC) and SOLV (water + ions) groups; τT = 1 ps for each |
| Pressure coupling | C-rescale, semiisotropic; 1 bar in both coupled directions; τP = 5 ps; compressibility 4.5 × 10⁻⁵ bar⁻¹ in each direction |
| Nonbonded list | Verlet; `rlist = 1.2 nm`; buffer tolerance 0.005 kJ mol⁻¹ ps⁻¹ |
| Van der Waals interactions | Force switch from 1.0 to 1.2 nm; dispersion correction off |
| Electrostatics | Particle-mesh Ewald; real-space cutoff 1.2 nm; Fourier spacing 0.12 nm; interpolation order 4; relative tolerance 10⁻⁵ |
| Constraints | Hydrogen-containing bonds with LINCS, order 4 and two iterations; rigid TIP3P with SETTLE |
| Production restraints | Position and dihedral restraint macros disabled |
| Center-of-mass motion removal | Whole-system linear motion removed every 100 steps |
| Coordinate output | Compressed coordinates every 5,000 steps (10 ps); precision 1,000 nm⁻¹ |
| Energy output | Every 500 steps (1 ps) |

The [production MDP](configs/simulation/production.mdp) is identical across the original six systems. Parallel execution and the neighbor-list runtime setting are described in [computational environments](ENVIRONMENT.md).

## Directed V368R–ascr3 starting poses

Both directed starts use the equilibrated `npt05` receptor and membrane of `v368r_ascr3_lid_rebuilt_own_pose_v1`, retaining its rebuilt lid. The Arg368-directed design changes the ligand pose and Arg368 side chain; the Arg189-directed design changes only the ligand pose. These are two alternative initial conditions, each followed by one 100 ns trajectory.

Starting-pose selection records and designed/minimized structures are external data. A new calculation must supply compatible structures and retain candidate identity, atom mapping, donor–H–acceptor geometry, internal strain and assembly records. A designed pre-minimization geometry must not be reported as a production-trajectory observation.

The system configuration specifies 6,408 receptor atoms plus 46 ligand atoms, with 145,881 atoms in the full membrane system. Directed starts use fresh velocities and the same heating, hold and five-stage NPT schedule. Their [MDP files](configs/simulation/polar-starts/) and [system mapping](configs/simulation/systems.json) specify seed and preparation settings separately. Input preparation and the resulting atom order must match these settings. The production integration and physical parameters in the table above apply to these systems; runtime settings are recorded in [ENVIRONMENT.md](ENVIRONMENT.md).

## Inputs, atom identity and processing

Each system contains one production trajectory spanning 0–100,000 ps at 10 ps intervals: 10,001 unique timepoints. Repeated timestamps at 1,030 ps for N84V–ascr18 and 26,970 ps for rebuilt-lid WT–ascr3 are deduplicated only when their coordinates agree. Statistical denominators use unique timepoints.

Each system used its own full topology, with 406 receptor residues and a separate ligand molecule block. Atom correspondence used residue 84/368 identities, atom names, residue numbers and atom order. Reduction preserved receptor and ligand order. GROMACS PDB output truncates ASCR18/ASCR3 to ASCR; ligand selection therefore uses the fixed atom indices in `configs/systems.json`. The `partner.selection`, `partner.com_selection` and `selection_rmsd_selection` fields define the ligand, its COM heavy atoms and its RMSD heavy atoms, respectively. Atom indices are zero-based; residue numbers are one-based original receptor numbers.

Processing used GROMACS 2025.4: concatenate full-system segments, `-pbc mol -center` around the receptor while retaining receptor and ligand, `-pbc nojump`, then receptor-backbone `-fit rot+trans`. Ligand motion relative to the receptor was preserved. Coordinate diagnostics covered finite coordinates, the complete time grid, heavy-atom covalent lengths (3.0 Å screen), neighboring receptor CA distances (4.5 Å) and 10 ps relative ligand COM steps (10 Å).

## Basic metrics

The primary descriptive window is 50–100 ns inclusive. Full curves, the initial portion and nonoverlapping 10 ns blocks describe the time dependence within each single trajectory. The last block includes the 100 ns endpoint. Means, occupancies and block summaries describe repeated samples within each trajectory.

| Metric | Definition | Sampling / primary denominator |
|---|---|---|
| Backbone RMSD, CA Rg | Receptor backbone fitted to its own production frame 0; receptor CA radius of gyration | 10 ps; 5,001 late frames |
| Ligand RMSD | All ligand heavy atoms after equal-weight receptor backbone Kabsch fit to its own frame 0 | 10 ps; 5,001 late frames |
| Geometric H-bond | N/O/S donor candidates, frame-0 nearest-H pairing within 1.35 Å, donor–acceptor ≤3.0 Å and D–H–A ≥150° | 500 ps; 101 late samples |
| COM distances | Mass-weighted ligand heavy-atom COM to whole receptor or its own initial 4 Å residue cluster | 10 ps; 5,001 late frames |
| CA RMSF | Receptor CA fluctuations about the 50–100 ns mean on backbone-fitted processed coordinates | 5,001 late frames; 406 CA positions |
| Buried surface area | `SASA(receptor) + SASA(ligand) − SASA(complex)`, Shrake–Rupley with a 1.4 Å probe and 240 sphere points; no division by two | 500 ps; 201 full / 101 late samples |
| Initial contact retention | Minimum receptor-residue / ligand-residue heavy-atom distance ≤4 Å, divided by the number of own frame-0 contacts | 500 ps; 201 full / 101 late samples |

Distances and RMSD/RMSF are in Å; BSA is Å²; XTC time is ps and reported time is ns. MDTraj SASA nm² is multiplied by 100. Each system has a distinct initial contact cluster and heavy-atom ligand set. Normalized BSA uses its own initial value.

The common-frame ligand paths additionally fit each receptor backbone (`protein and backbone`) to WT–ascr18 frame 0, then use heavy-atom geometric ligand centroids and all-residue 4 Å contacts. Sampling is 500 ps, with 100 early samples below 50 ns and 101 late samples from 50–100 ns inclusive. Full sampled paths and each arm's own initial ligand pose are retained. This common coordinate frame displays the eight prepared starting poses and their subsequent motion.

## Local fits and selection definitions

Local analysis samples every 100 ps: 1,001 full timepoints, 501 in 50–100 ns. Only frame-0 data define fits. The local parameters are in [mechanism.json](configs/mechanism.json), with `stride_frames = 10` and `window_ns = [50, 100]`. Fixed selections are in [selections.json](configs/selections.json):

- In each full initial membrane system, POPC P atoms are split by median z into two leaflets; mean leaflet z coordinates delimit an operational membrane slab. Receptor frame-0 CA positions within this slab are intersected with simplified DSSP H residues. The intersection across six systems, retaining contiguous runs of at least seven residues, gives 123 common core residues. Simplified H includes alpha, 3₁₀ and pi helices.
- Pocket residues are the union across all six initial frames of receptor residues within 6 Å of ligand heavy atoms, excluding lid 164–178: 37 common residues.
- Global, core and pocket fits use N/CA/C/O atoms and each system's own frame 0. Outputs retain ligand heavy-atom RMSD, ligand geometric-centroid displacement and lid-backbone RMSD for all three fits. Lid–pocket distance measures the separation between their backbone geometric centroids.

`configs/selections.json` retains the fixed analysis residue sets and their derivation definition. Per-system frame-0 observations are external data. These same 123 core and 37 pocket residues define the comparison for the two directed starts, with each trajectory fitted to its own frame 0. Users can derive new selections with `--prepare`; changed selections define a different comparison.

## Local chemical and spatial descriptors

For residues 84 and 368, side-chain heavy atoms are compared to all other residue heavy atoms at cutoffs 3.5, 4.0 and 4.5 Å. Sequence separation is retained. Chi1 uses N–CA–CB–CG (Val uses CG1); circular distance assigns residue-specific rotamer bins centered at −60°, +60° and 180° and handles the ±180° periodic boundary.

Residues 84, 368, 102, 186, 189 and 345 have atom-resolved ligand minimum distances, nearest atom names and ≤4.5 Å atom pairs. Ligand groups derive from each ligand's topology bonds: carboxylate carbon and its singly bonded oxygens, the oxygen-containing six-membered sugar ring and adjacent oxygens / terminal methyl, and remaining side chain. Each `ligand_atom_mapping.csv` gives that ligand's names, types, charges, indices and bond neighbors.

Local candidate hydrogen bonds use topology donor–H bonds. Donors are H-bearing N/O/S; acceptors are protein O atoms, unprotonated histidine ring N and ligand O atoms. Geometry remains donor–acceptor ≤3 Å and D–H–A ≥150°. Per-residue occupancy counts the union of qualifying direct receptor–ligand events at each timepoint once. This local chemical definition differs from the coarse basic driver.

Two-step spatial adjacency connects site side chain → intermediate residue → anchor, both edges ≤4.5 Å, with sequence gaps greater than two for both edges. Joint occupancy measures simultaneous geometric contacts and additionally requires anchor–ligand ≤4 Å at the same time. The first sustained event is the start of the first qualifying interval spanning at least 1 ns (11 samples at 100 ps). Events already present at t=0 are identified as initial contacts.

Atom-pair analysis measures receptor contacts for 84–25, 84–28, 84–80, 84–90, 84–171, 84–186, 84–189, 368–38, 368–70, 368–77, 368–364, 368–345 and 364–345, together with Arg368 chi2/3/4. The charged-end Arg368–Asp70 distance is NE/NH1/NH2 to OD1/OD2, distinct from whole-side-chain minimum distance. Contact events use distances ≤4 Å; separation events use distances >4 Å. Sustained-event times use the unsmoothed 100 ps samples. The timeline figure displays a 1 ns moving average. Basic and local occupancies are compared on the same 500 ps sampling grid.

References: [GROMACS 2025.4 simulation parameters](https://manual.gromacs.org/2025.4/user-guide/mdp-options.html); [GROMACS trajectory tools](https://manual.gromacs.org/2025.4/onlinehelp/gmx-trjconv.html); [MDAnalysis alignment](https://docs.mdanalysis.org/2.10.0/documentation_pages/analysis/align.html); [MDTraj Shrake–Rupley SASA](https://mdtraj.readthedocs.io/en/stable/api/generated/mdtraj.shrake_rupley.html); [MDTraj DSSP](https://mdtraj.readthedocs.io/en/stable/api/generated/mdtraj.compute_dssp.html).

## Polar contacts and transition geometry

Targeted polar analysis reads all 10,001 frames at 10 ps spacing in the five systems listed in [polar-anchors.json](configs/polar-anchors.json). It measures Arg368, Arg189 and Arg345 against the ligand carboxylate and the competing Arg368–Asp70 contact. Salt bridges use minimum charged N–O distance ≤4.0 Å, with 3.5/4.5 Å sensitivity thresholds. Short linear hydrogen bonds require donor–acceptor ≤3.2 Å, H–acceptor ≤2.4 Å and D–H–A ≥150° simultaneously, using actual topology donor–H bonds. This polar definition has its own outputs and differs from the 3.0 Å basic/local definition above. A union of eligible atom pairs counts each frame once; chemical absence remains missing rather than measured zero.

The [all-residue screen](configs/polar-screen.json) covers all 406 receptor residues in both directed starts, WT–ascr3 and V368R–ascr18. Features distinguish receptor residue, side chain/backbone, ligand chemical group and donor direction. It includes both N/O hydrogen-bond directions, Arg/Lys/charged His and protonated N-terminal cations against ligand carboxylate. Ligand groups follow bonds rather than shared atom labels. [Screen parameters](configs/polar-screen-definition.json) specify primary, tight (3.0 Å / 2.4 Å / 150°) and relaxed (3.5 Å / 2.6 Å / 135°) hydrogen bonds. These are direct contacts in receptor–ligand coordinates; water-mediated interactions and interaction energies are outside these definitions.

The primary window is 50–100 ns inclusive (5,001 frames). A pattern requires target occupancy ≥70%, occupancy ≤20% in each control, and at least two adjacent 10 ns target blocks at ≥70%. Additional 20 ns windows begin every 10 ns from 0 through 80 ns; only windows ending at 100 ns include their right endpoint. Inter-residue combinations are simultaneous AND, with exact frequency pruning and no order cap; a target pattern must contain a non-Arg345 residue. Full single-feature and combination tables retain failures, controls and all thresholds. Uninterrupted spans do not fill gaps; a separate duration permits gaps up to 0.05 ns. Overlapping 20–40 and 30–50 ns observations describe one transient stage, not independent confirmations.

[Asp89 transition parameters](configs/polar-transition.json) define the O3–bonded H···Asp89 OD1/OD2 geometry and phases 15–18, 20–25, 45–48, 51–54 and 60–80 ns, each right-exclusive. These windows describe an already observed event. Core fitting uses the fixed common core; local fitting uses residues 82–95. Every 100 ps, all cross-window frame pairs contribute a symmetric two-endpoint decomposition of the O3–nearest Asp89 oxygen distance change. Averaging the two movement orders assigns donor and acceptor coordinate contributions whose sum equals the total distance change. This describes motion, not forces or causality. Asp89 χ2 is stored using OD1; the equivalent OD1/OD2 oxygens require 180° periodicity when interpreting rotamers.

The transition analysis also retains side-chain dihedrals, Arg345/Arg189 ligand distances, Arg368–Asp70 geometry and all-residue neighbor distances every 100 ps. A one-intermediate proximity requires both site-to-residue edges ≤4.5 Å simultaneously, excluding ±2 sequence neighbors of each endpoint. Actual phase structures are majority-state frames nearest the local-fit mean sugar/Asp oxygen coordinates. Fixed compatibility frames are the nearest hydrogen-bonded frames to 22.5, 27.5, 32.5, 37.5 and 43.38 ns within 0.5 ns. The selected frame indices, times and file hashes accompany the structures.

## Static ligand compatibility

[Compatibility parameters](configs/pose-compatibility.json), [chemical graph and stereochemistry](configs/ligand_geometry_contract.json), [ascr3](configs/ascr3.sdf) and [ascr18](configs/ascr18.sdf) define the comparison in five actual Arg189-directed receptor frames. The receptor, shared sugar heavy atoms and attached hydrogens remain fixed; remaining ligand atoms can move. The result is a constructed static model, not an observed ascr18 trajectory frame.

RDKit ETKDGv3 uses seed 20261016, at most 96 conformers per ligand and one thread. After MMFF94 internal optimization, the 16 least-overlapping fixed-scaffold candidates are refined for each pocket; ascr3 additionally includes the supplied trajectory pose. L-BFGS-B runs at most 700 iterations with MMFF94 internal energy and a steric penalty of 100 times the squared overlap beyond the summed heavy-atom radii minus 0.4 Å. Optimization uses protein atoms within 26 Å of the scaffold center; final geometry uses the whole protein.

A passing candidate must converge, have no protein–ligand heavy-atom pair below 2.2 Å, maximum van der Waals overlap ≤0.8 Å, maximum bond-length error ≤0.18 Å, unchanged stereochemistry/double-bond geometry and internal strain ≤12 kcal/mol relative to the same ligand's lowest free-chain energy with that fixed scaffold. Full candidate tables retain failures and optimizer status. Selected PDB and chemically explicit SDF files preserve all five pockets for both ligands, including selections from pockets with no passing candidate. Solvent, membrane and receptor relaxation are absent; these scores describe fixed-pocket geometric compatibility and cannot rank binding affinity or predict dynamic persistence.
