import re, os, json, argparse

import numpy as np
import seaborn as sns
from ase.io import read
import matplotlib.pyplot as plt


def get_e_from_lammps_log(lammps_log):
    f = open(lammps_log, "r")
    Lines = f.readlines()
    patten = r"(\d+\s+\-+\d*\.?\d+)"
    e_pot = []
    for i, line in enumerate(Lines):
        s = line.strip()
        match = re.match(patten, s)
        if match != None:
            D = np.fromstring(s, sep=" ")
            e_pot.append(D[1])
    return np.array(e_pot)


def plot_energy_trajectories(
    energy_list, delta=False, label_energy=None, save=False, save_folder=None
):
    """
    Given a list of lists of energies, plot them
    """
    assert (
        delta == False or label_energy == None
    ), "Cannot plot delta energy with label_energy"

    for energy in energy_list:
        if delta:
            energy = energy - energy[0]
        plt.plot(energy)

    plt.xlabel("Opt step")
    plt.title("Energy vs. opt step (n={})".format(len(energy_list)))

    if delta:
        plt.ylabel("$\Delta$ Energy (eV)")
    else:
        plt.ylabel("Energy (eV)")

    if label_energy != None:
        plt.axhline(y=label_energy, color="red", label="DFT init energy")
        plt.legend()

    plt.show()
    if save:
        final_name = "final_energy" if not delta else "final_energy_delta"
        file_save = os.path.join(save_folder, "{}.png".format(final_name))
        plt.savefig(file_save, dpi=300)


def plot_hist_final_energy(
    energy_list, label_energy=None, save=False, save_folder=None
):
    """
    Given a list of energy trajectories, plot the histogram of final energies
    """
    energies_final = []
    for energy in energy_list:
        energies_final.append(energy[-1])
    # make kde plot
    # sns.kdeplot(energies_final, fill=True, bw_adjust=0.3)
    sns.histplot(energies_final, kde=True, stat="density", bins=30)
    plt.xlabel("Final energy (eV)")
    plt.ylabel("Density")
    plt.title("Final energy distribution (n={})".format(len(energy_list)))
    if label_energy != None:
        plt.axvline(x=label_energy, color="red", label="DFT init energy")
        plt.legend()

        lower_count = 0
        for energy in energy_list:
            final_e = energy[-1]

            if final_e < float(label_energy):
                lower_count += 1
        print(
            "% trajectories w/ final E < than DFT init E: {:.2f}".format(
                lower_count / len(energy_list) * 100
            )
        )

    # plt.hist(energies_final)
    plt.show()
    if save:
        file_save = os.path.join(save_folder, "final_energy_hist.png")
        plt.savefig(file_save, dpi=300)


def main():
    # argparse to options json
    parser = argparse.ArgumentParser()
    parser.add_argument("--options", type=str, default="./viz_options.json")
    args = parser.parse_args()
    with open(args.options) as f:
        options = json.load(f)

    root = options["root"]
    folder_save_best_structs = options["folder_save_best_structs"]
    save_best_structs = options["save_best_structs"]
    save_lowest_per_trajectory = options["save_last_in_ea_md"]
    label_energy = options["label_energy"]

    # get all subfolders
    subfolders = [f.path for f in os.scandir(root) if f.is_dir()]
    print("found {} subfolders".format(len(subfolders)))
    # get all trajectories

    energy_list = []
    flat_energy_list = []
    flat_force_list = []
    index_list = []

    for folder_ind, subfolder in enumerate(subfolders):
        # open optimized.traj and opt.traj in each subfolder
        opt_traj = os.path.join(subfolder, "opt.traj")  # trajectory
        optimized_traj = os.path.join(subfolder, "optimized.traj")  # final
        #lammps_traj = os.path.join(subfolder, "md.lammpstrj")
        #lammps_log = os.path.join(subfolder, "log.lammps")

        # check that all files exist in subfolder
        if (
            os.path.exists(opt_traj)
            #and os.path.exists(lammps_traj)
            #and os.path.exists(lammps_log)
            and os.path.exists(optimized_traj)
        ):
            # read in the optimized.traj
            atoms_optimized = read(optimized_traj, index=":")
            # read in the opt.traj
            atoms_opt = read(opt_traj, index=":")
            # read md.lammpstrj
            # images = read(lammps_traj, ":")  # doesn't contain forces

            e_pot = []
            # get energy method 1
            for index, atom in enumerate(atoms_opt):
                energy = atom.get_potential_energy()
                forces = atom.get_forces()
                # print(f"Energy of trajectory {index+1}: {energy} eV")
                e_pot.append(energy)

                if save_best_structs:
                    flat_energy_list.append(energy)
                    flat_force_list.append(forces)
                    index_list.append((folder_ind, index))

            e_pot = np.array(e_pot)

            # get energy method 2
            #e_pot = get_e_from_lammps_log(lammps_log)

            energy_list.append(e_pot)
            # print("num energies: {}".format(len(e_pot)))

    plot_energy_trajectories(energy_list, delta=True, save=True, save_folder=root)
    plot_energy_trajectories(
        energy_list, delta=False, label_energy=label_energy, save=True, save_folder=root
    )
    plot_hist_final_energy(
        energy_list, label_energy=label_energy, save=True, save_folder=root
    )

    if save_best_structs:
        # get top n lowest energy structures
        n = 100
        # get indices of lowest n energies
        lowest_n_indices = np.argsort(flat_energy_list)[:n]
        # get corresponding indices of subfolders and images
        lowest_n_full_indices = [index_list[i] for i in lowest_n_indices]

        for ind, target_index in enumerate(lowest_n_full_indices):
            print(
                "folder ind {} frame {} energy {}".format(
                    target_index[0],
                    target_index[1],
                    flat_energy_list[lowest_n_indices[ind]],
                )
            )
            # if folder_save_best_structs doesn't exist, create it
            if not os.path.exists(folder_save_best_structs):
                os.makedirs(folder_save_best_structs)
            target_folder = subfolders[target_index[0]]
            # print(target_folder)
            opt_traj = os.path.join(target_folder, "opt.traj")
            # write to new folder
            write_path_traj = os.path.join(
                folder_save_best_structs, "best-{}.traj".format(ind)
            )
            write_path_xyz = os.path.join(
                folder_save_best_structs, "best-{}.xyz".format(ind)
            )
            atoms_opt = read(opt_traj, index=":")
            atoms_opt[target_index[1]].write(write_path_traj)
            atoms_opt[target_index[1]].write(write_path_xyz)
            # add energies to target_index
            target_index += (flat_energy_list[lowest_n_indices[ind]],)
            lowest_n_full_indices[ind] = target_index
            # save dictionary of lowest_full_indices
            with open(
                os.path.join(folder_save_best_structs, "lowest_n_indices.json"), "w"
            ) as f:
                json.dump(lowest_n_full_indices, f)

    if save_lowest_per_trajectory:
        for folder_ind, subfolder in enumerate(subfolders):
            # print(subfolder)
            # open optimized.traj and opt.traj in each subfolder
            opt_traj = os.path.join(subfolder, "opt.traj")  # trajectory
            optimized_traj = os.path.join(subfolder, "optimized.traj")  # final
            #lammps_traj = os.path.join(subfolder, "md.lammpstrj")
            #lammps_log = os.path.join(subfolder, "log.lammps")

            # check that all files exist in subfolder
            if (
                os.path.exists(opt_traj)
                #and os.path.exists(lammps_traj)
                #and os.path.exists(lammps_log)
                and os.path.exists(optimized_traj)
            ):
                atoms_opt = read(opt_traj, index="-1")
                print("final opt energy: {}".format(atoms_opt.get_potential_energy()))
                folder_number = int(subfolder.split("_")[-1])
                write_path_traj = os.path.join(
                    folder_save_best_structs, "last-{}.traj".format(folder_number)
                )
                write_path_xyz = os.path.join(
                    folder_save_best_structs, "last-{}.xyz".format(folder_number)
                )
                atoms_opt.write(write_path_traj)
                atoms_opt.write(write_path_xyz)


main()
