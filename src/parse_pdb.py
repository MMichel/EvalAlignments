import sys
import operator
import numpy as np
from collections import defaultdict


def parse_atm_record(line):

    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])
    
    return record


def read(pdbfile, chain='', model=1):

    header = ''
    res_lst = []
    atm_lst = []
    tail = ''

    seen_atoms = False
    in_atoms = False
    in_model = False
    curr_resi = 0
    prev_resi = 0
    
    for line in pdbfile:
        if line.startswith('MODEL') and int(line.strip().split()[-1]) == model:
            in_model = True
            header += line
        elif in_model and line.startswith('TER'):
            atm_record = parse_atm_record(atm_lst[-1])
            if chain and not chain == atm_record['chain']:
                continue
            seen_atoms = True
            #print "seen_atoms model"
            #print len(res_lst)
            #in_atoms = False
            tail += line
        elif in_model and line.startswith('ENDMDL'):
            in_model = False
            tail += line
        elif line.startswith('MODEL') and not int(line.strip().split()[-1]) == model:
            continue
        elif not in_model and line.startswith('TER'):
            atm_record = parse_atm_record(atm_lst[-1])
            if chain and not chain == atm_record['chain']:
                continue
            seen_atoms = True
            #in_atoms = False
            #print "seen atoms"
            #print len(res_lst)
            continue
        elif not in_model and line.startswith('ENDMDL'):
            continue
        elif not line.startswith('ATOM') and not seen_atoms:
            header += line
        elif not line.startswith('ATOM') and seen_atoms:
            tail += line
        elif in_model or ((not seen_atoms) and (not in_model)):
            atm_record = parse_atm_record(line)
            if chain and not chain == atm_record['chain']:
                atm_lst = [line]
                continue
            if not in_atoms:
                curr_resi = atm_record['res_no']
                prev_resi = curr_resi
            in_atoms = True
            curr_resi = atm_record['res_no']
            if curr_resi == prev_resi:
                atm_lst.append(line)
            else:
                res_lst.append(atm_lst)
                atm_lst = [line]
            prev_resi = curr_resi
    res_lst.append(atm_lst)
     
    pdbfile.close()
    pdb_lst = [header, res_lst, tail]
    return pdb_lst


def read_chain(pdbfile, chain):

    header = ''
    res_lst = []
    atm_lst = []
    tail = ''

    seen_atoms = False
    curr_resi = 0
    prev_resi = 0
    
    for line in pdbfile:
        if not line.startswith('ATOM') and not seen_atoms:
            header += line
        elif not line.startswith('ATOM') and seen_atoms:
            tail += line
        else:
            atm_record = parse_atm_record(line)
            if not atm_record['chain'] == chain:
                continue
            if not seen_atoms:
                curr_resi = atm_record['res_no']
                prev_resi = curr_resi
            seen_atoms = True
            curr_resi = atm_record['res_no']
            if curr_resi == prev_resi:
                atm_lst.append(line)
            else:
                #atm_lst.append(line)
                res_lst.append(atm_lst)
                atm_lst = [line]
            prev_resi = curr_resi
    res_lst.append(atm_lst)
     
    pdbfile.close()
    pdb_lst = [header, res_lst, tail]
    return pdb_lst



def write(pdb_lst, outfile):

    outfile.write(pdb_lst[0])

    for res in pdb_lst[1]:
        for atm in res:
            outfile.write(atm)
            
    outfile.write(pdb_lst[2])
    outfile.close()


def get_coordinates(pdbfile, chain):

    res_dict = defaultdict(list)

    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        if atm_record['chain'] != ' ' and atm_record['chain'] != chain:
            continue

        res_i = atm_record['res_no']
        atm = [atm_record['x'], atm_record['y'], atm_record['z']]

        res_dict[res_i].append(np.array(atm))
        
    pdbfile.close()
    return sorted(res_dict.iteritems(), key=operator.itemgetter(0))


def get_res_dict(pdbfile, chain):

    cb_lst = []
    res_dict = defaultdict(list)

    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue

        atm_record = parse_atm_record(line)

        if atm_record['chain'] != ' ' and atm_record['chain'] != chain:
            continue

        res_i = atm_record['res_no']
        
        if atm_record['insert'] == 'X':
            res_i = res_i * 0.001

        atm = [float('inf'), float('inf'), float('inf')]

        if atm_record['atm_name'] == 'CA':
                atm = [atm_record['x'], atm_record['y'], atm_record['z']]
                res_dict[res_i].append(np.array(atm))   
        elif atm_record['atm_name'] == 'CB':
                atm = [atm_record['x'], atm_record['y'], atm_record['z']]
                res_dict[res_i].append(np.array(atm)) 
    
    return res_dict


def get_ca_coordinates(pdbfile, chain):

    res_dict = get_res_dict(pdbfile, chain)

    ca_lst = []

    # need to sort to get the sequence correct
    sorted_keys = sorted(res_dict.keys())
    
    for i in sorted_keys:
        ca_lst.append(res_dict[i][0])
    pdbfile.close()
    return ca_lst


def get_cb_coordinates(pdbfile, chain):

    res_dict = get_res_dict(pdbfile, chain)

    cb_lst = []
    tmp_i = 0

    # need to sort to get the sequence correct
    sorted_keys = sorted(res_dict.keys())
    
    for i in sorted_keys:
        if len(res_dict[i]) > 1:
            tmp_i += 1
            cb_lst.append(res_dict[i][-1])
        elif len(res_dict[i]) == 1:
            tmp_i += 1
            cb_lst.append(res_dict[i][0])
    #print atm_count 
    pdbfile.close()
    return cb_lst


def get_atom_seq(pdbfile, chain='', model=1):

    three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'UNK': 'X'}
    res_dict = {}
    
    in_model = False
 
    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    res_name = ''
    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        if atm_record['chain'] != ' ' and atm_record['chain'] != chain:
            continue
        if atm_record['atm_name'] != 'CA':
            continue

        res_i = atm_record['res_no']
        if atm_record['insert'] == 'X':
            res_i = res_i * 0.001
         
        #print atm_record['res_name']
        if atm_record['res_name'] in three_to_one:
            #res_name = three_to_one[atm_record['res_name']]
            #print res_name
            res_name = three_to_one[atm_record['res_name']]
        #else:
            #res_name = ''
            #continue

        res_dict[res_i] = res_name

    res_lst = sorted(res_dict.iteritems(), key=operator.itemgetter(0))
    atom_seq = ''

    for res in res_lst:
        atom_seq += res[1]

    pdbfile.close()
    return atom_seq


def get_first_chain(pdbfile):

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        break

    return atm_record['chain']
 

def get_acc(pdbfile):

    for line in pdbfile:
        if line.startswith('HEADER'):
            return line[62:66].lower()
    # if no header line in pdb file:
    return ''


if __name__ == '__main__':

    pdbfile = open(sys.argv[1], 'r')
    chain = sys.argv[2]
    #print get_atom_seq(pdbfile, chain)
    pdbfile.close()
    #pdbfile = open(sys.argv[1], 'r')
    #print get_coordinates(pdbfile)
    #pdbfile.close()
