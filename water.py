import os
import pathlib
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename

pd.set_option('display.max_columns', None)


def process_hb2():
    content = []

    Tk().withdraw()
    hb2_name = askopenfilename()
    # Check if the file is a .hb2 file
    if pathlib.Path(hb2_name).suffix == ".hb2":
        with open(hb2_name, 'r') as hb2:
            # fill each line as an element of the list "content"
            content = hb2.readlines()
            # Remove \n from the end of each line in the list
            content = [line.strip() for line in content]
    else:
        print("The file is not in the right format (.hb2 file)")

    # Skip comments (line 0 to 7)
    content = content[8:]

    # An example of an entry in .hb2 file, read top down for the column number in
    # Table 1 of https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html
    #
    # 00000000011111111112222222222333333333344444444445555555555666666666677777777778
    # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
    # A0069-VAL N   A2001-HOH O   3.27 MH  -2 -1.00 163.8  2.28  -1.0  -1.0     1
    #
    # donor           :  A0069 (whitespace removed)
    # donor_amino     :  VAL   (whitespace removed)
    # donor_type      :  N     (whitespace removed)
    # acceptor        :  A2001 (whitespace removed)
    # acceptor_amino  :  HOH   (whitespace removed)
    # acceptor_type   :  O     (whitespace removed)
    # da_dist         :  3.27
    # cat_da          :  MH    (whitespace removed)
    # ... so on, refer to Table 1 ...

    donor, donor_amino, donor_type, acceptor, acceptor_amino, acceptor_type,\
        da_dist, cat_da, num_aas, dist, dha_angle, ha_dist, h_a_aa_angle, d_a_aa_angle = ([] for _ in range(14))

    for line in content:
        # Since python starts indexing at zero, the start index is always the start position -1
        # Since the last index is not included in a range of indices, we don't need to subtract 1 from the end position
        donor.append(line[0:5])
        donor_amino.append(line[6:9].strip())  # Remove white spaces
        donor_type.append(line[9:13].strip())  # Remove white spaces
        acceptor.append(line[14:19])
        acceptor_amino.append(line[20:23].strip())  # Remove white spaces
        acceptor_type.append(line[23:27].strip())  # Remove white spaces
        da_dist.append(float(line[27:32]))  # Change string to float
        cat_da.append(line[33:35])
        num_aas.append(float(line[36:39]))  # Change string to float
        dist.append(float(line[40:45]))  # Change string to float
        dha_angle.append(float(line[46:51]))  # Change string to float
        ha_dist.append(float(line[52:57]))  # Change string to float
        h_a_aa_angle.append(float(line[58:63]))  # Change string to float
        d_a_aa_angle.append(float(line[64:69]))  # Change string to float

    # Create a dictionary to be changed to a dataframe
    keys = ['donor', 'donor_amino', 'donor_type', 'acceptor', 'acceptor_amino', 'acceptor_type',
            'da_dist', 'cat_da', 'num_aas', 'dist', 'dha_angle', 'ha_dist', 'h_a_aa_angle', 'd_a_aa_angle']
    values = [donor, donor_amino, donor_type, acceptor, acceptor_amino, acceptor_type,
              da_dist, cat_da, num_aas, dist, dha_angle, ha_dist, h_a_aa_angle, d_a_aa_angle]
    hb2_dict = dict(zip(keys, values))
    hb2_df = pd.DataFrame(hb2_dict)
    # print(hb2_df[:5])

    return hb2_df


def is_nucleotide(x):
    return x in ['A', 'C', 'T', 'U', 'G', 'ATP']  # Standard Nucleotides: C, A, U, G, T, also ATP


def is_amino(x):
    return x in ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE',
                 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


def find_water_bridges(hb2_df, double=False):
    sol = hb2_df.loc[(result['donor_amino'] == 'SOL') | (hb2_df['acceptor_amino'] == 'SOL')]
    # print(sol.head())
    # print(len(sol))

    sol_nuc_acc, sol_amino_acc, sol_nuc_don, sol_amino_don, sol_sol = ([] for _ in range(5))

    for i in range(len(sol)):
        if 'SOL' in sol['donor_amino'].iloc[i]:
            if is_nucleotide(sol['acceptor_amino'].iloc[i]):
                if [sol['donor'].iloc[i], sol['acceptor'].iloc[i]] not in sol_nuc_acc:
                    sol_nuc_acc.append([sol['donor'].iloc[i], sol['acceptor'].iloc[i]])
            elif is_amino(sol['acceptor_amino'].iloc[i]):
                if [sol['donor'].iloc[i], sol['acceptor'].iloc[i]] not in sol_amino_acc:
                    sol_amino_acc.append([sol['donor'].iloc[i], sol['acceptor'].iloc[i]])

        elif 'SOL' in sol['acceptor_amino'].iloc[i]:
            if is_nucleotide(sol['donor_amino'].iloc[i]):
                if [sol['donor'].iloc[i], sol['acceptor'].iloc[i]] not in sol_nuc_don:
                    sol_nuc_don.append([sol['donor'].iloc[i], sol['acceptor'].iloc[i]])
            elif is_amino(sol['donor_amino'].iloc[i]):
                if [sol['donor'].iloc[i], sol['acceptor'].iloc[i]] not in sol_amino_don:
                    sol_amino_don.append([sol['donor'].iloc[i], sol['acceptor'].iloc[i]])

    # print(sol_nuc_acc, '\n')
    # print(sol_amino_don, '\n')
    # print(sol_nuc_don, '\n')
    # print(sol_amino_acc, '\n')

    swb = []
    for pair_nuc in sol_nuc_don:
        sol_id = pair_nuc[1]

        for pair_aa in sol_amino_acc:
            if pair_aa[0] == sol_id:
                swb.append([pair_nuc[0], sol_id, pair_aa[1], "don - acc"])

        for pair_amino in sol_amino_don:
            if pair_amino[1] == sol_id:
                swb.append([pair_nuc[0], sol_id, pair_amino[0], "don - don"])

    for pair_nuc in sol_nuc_acc:
        sol_id = pair_nuc[0]

        for pair_aa in sol_amino_don:
            if pair_aa[1] == sol_id:
                swb.append([pair_nuc[1], sol_id, pair_aa[0], "acc - don"])

        for pair_amino in sol_amino_acc:
            if pair_amino[0] == sol_id:
                swb.append([pair_nuc[1], sol_id, pair_amino[1], "acc - acc"])

    swb_df = pd.DataFrame(swb, columns=['Nucleo', 'Water', 'A.acid', 'Bonding Type'])

    if double is True:
        for i in range(len(sol)):
            if 'SOL' in sol['donor_amino'].iloc[i] and 'SOL' in sol['acceptor_amino'].iloc[i]:
                if [sol['donor'].iloc[i], sol['acceptor'].iloc[i]] not in sol_sol:
                    sol_sol.append([sol['donor'].iloc[i], sol['acceptor'].iloc[i]])

        # print(sol_sol)

        dwb = []
        for pair_sol in sol_sol:
            sol_id_don, sol_id_acc = pair_sol[0], pair_sol[1]

            for pair_nuc in sol_nuc_don:
                if pair_nuc[1] == sol_id_don:
                    for pair_aa in sol_amino_acc:
                        if pair_aa[0] == sol_id_acc:
                            dwb.append([pair_nuc[0], sol_id_don, sol_id_acc, pair_aa[1], "don - acc"])

                    for pair_amino in sol_amino_don:
                        if pair_amino[1] == sol_id_acc:
                            dwb.append([pair_nuc[0], sol_id_don, sol_id_acc, pair_amino[0], "don - don"])

            for pair_nuc in sol_nuc_acc:
                if pair_nuc[0] == sol_id_acc:
                    for pair_aa in sol_amino_don:
                        if pair_aa[1] == sol_id_don:
                            dwb.append([pair_nuc[1], sol_id_acc, sol_id_don, pair_aa[0], "acc - don"])

                    for pair_amino in sol_amino_acc:
                        if pair_amino[0] == sol_id_don:
                            dwb.append([pair_nuc[1], sol_id_acc, sol_id_don, pair_amino[1], "acc - acc"])

        dwb_df = pd.DataFrame(dwb, columns=['Nucleo', 'Water', 'Water', 'A.acid', 'Bonding Type'])
        dwb_df = dwb_df.shift()[1:]

        return swb_df, dwb_df

    return swb_df

result = process_hb2()
# print(result.head())
# print(len(result))

s_water, d_water = find_water_bridges(result, double=True)

print(s_water)
print(d_water)
