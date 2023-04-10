# Ruize Liu
# GRS analysis 
# Calcuate GRS for each sample with PLINK


import os, argparse, sys

parser = argparse.ArgumentParser(description='plink --clump & --score function\n'
                                             'result: out/baitA/baitA_int.clump out/baitA/baitA_non.clump\n'
                                             '        out/baitA/baitA_int.profile out/baitA/baitA_non.profile')
parser.add_argument('--plink', help='addr to plink e.g /home/user/bin/plink', default='plink')
parser.add_argument('--bfile', help='plink --bfile', action="store")
parser.add_argument('--clump', help='plink --clump', action="store")
parser.add_argument('--clump-kb', help='plink --clump-kb', action="store")
parser.add_argument('--clump-p1', help='plink --clump-p1', action="store")
parser.add_argument('--clump-p2', help='plink --clump-p2', action="store")
parser.add_argument('--clump-r2', help='plink --clump-r2', action="store")
parser.add_argument('--geno', help='plink --geno', action="store")
parser.add_argument('--keep', help='plink --keep', action="store")
parser.add_argument('--out', help='plink --out', action="store")
parser.add_argument('--maf', help='plink --maf', action="store")
parser.add_argument('--score', help='plink --score', action="store", nargs='+')
parser.add_argument('--tab', help='A table of file used for PRS', action="store")
parser.add_argument('--memory', help='plink --memory in MB', action="store")


# for LISA
def plink_fun(plink, bfile, clump, clump_kb, clump_p1, clump_p2, clump_r2, loc_tab, out, keep, geno, maf, score, memory):
    keep_cmd = ''
    geno_cmd = ''
    maf_cmd = ''
    memory_cmd = ''
    if keep:
        keep_cmd = '--keep {}'.format(keep)
    if geno:
        geno_cmd = '--geno {}'.format(geno)
    if maf:
        maf_cmd = '--maf {}'.format(maf)
    if memory:
        memory_cmd = '--memory {}'.format(memory)
    # clump
    ans = os.system(
        '{} --bfile {} --clump {} --clump-kb {} --clump-p1 {} --clump-p2 {} --clump-r2 {} --extract range {} --out {} {} {} {} {} --allow-no-sex'.format(
        plink, bfile, clump, clump_kb, clump_p1, clump_p2, clump_r2, loc_tab, out, keep_cmd, geno_cmd, maf_cmd, memory_cmd
        )
    )
    if ans != 0:
        print ans
        sys.exit()
    os.system('head -n 40 {}.log > {}.clump.log'.format(out, out))
    os.system('echo "...." >> {}.clump.log'.format(out))
    os.system('tail -n 20 {}.log >> {}.clump.log'.format(out, out))
    # prs
    ans = os.system(
        '{} --bfile {} --score {} --extract {} --out {} {} --allow-no-sex'.format(
            plink, bfile, score, out+'.clumped', out, memory_cmd
        )
    )
    if ans != 0:
        print ans
        sys.exit()
    os.system('mv {}.log {}.profile.log'.format(out, out))
    os.system('rm {}.nopred'.format(out))


args = parser.parse_args()
in_addr = args.tab
in_dir = os.path.split(in_addr)[0]


for line in open(in_addr):
    if line.startswith('#'):
        print line
        continue
    line_l = line.split()
    int_name = line_l[0]
    int_addr = os.path.join(in_dir, int_name + '_int.txt')
    non_addr = os.path.join(in_dir, int_name + '_non.txt')

    out_int_prefix = os.path.join(args.out, int_name + '_int')
    out_non_prefix = os.path.join(args.out, int_name + '_non')

    plink_fun(args.plink, args.bfile, args.clump, args.clump_kb, args.clump_p1, args.clump_p2, args.clump_r2,
              int_addr, out_int_prefix, args.keep, args.geno, args.maf, ' '.join(args.score), args.memory)

    plink_fun(args.plink, args.bfile, args.clump, args.clump_kb, args.clump_p1, args.clump_p2, args.clump_r2,
              non_addr, out_non_prefix, args.keep, args.geno, args.maf, ' '.join(args.score), args.memory)
