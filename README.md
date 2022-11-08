# Some python script for ab initio output

## imag\_freq
for remove imag freq

    option:
        -f: ab initio software
            gau: gaussian (default)
            bdf: BDF

        -s: scale add norm mode to coordinate (float)
            default: 0.1

    Result:
        2 * n xyz struct, initial coordinate plus scale norm mode 
        and minius scale norm mode. Where n is number of imag mode

    Usage:
        imag_freq.py xxx.log

    output:
        xxx00.xyz xxx01.xyz
