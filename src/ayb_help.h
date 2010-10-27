/**
 * \file ayb_help.h
 * Help text for AYB.
 */

//"12345678901234567890123456789012345678901234567890123456789012345678901234567890\n"
"\n"
"Options:\n"
"  -b  --blockstring <Rn[InCn...]>\n"
"\tHow to group cycle data in intensity files for analysis [required]\n"
"\t(R=Read, I=Ignore, C=Concatenate onto previous block)\n"
"\t(First R must precede first C)\n\n"
"  -x  --prefix <path>\t\tFilter input files by match to prefix [required]\n"
"  -c  --composition <prop. CG>\tGenome CG composition [default: 0.5]\n"
"  -d  --dataformat <format>\tInput format [default: txt]\n"
"\t\t\t\t(txt/cif)\n"
"  -f  --format <format>\t\tOutput format [default: fasta]\n"
"\t\t\t\t(fasta/fastq)\n"
"  -n  --niter <num>\t\tNumber of model iterations [default: 5]\n"
"  -m  --mu <num>\t\tAdjust range of quality scores [default: 1.0E-5]\n"
"\t\t\t\t(Smaller value for higher maximum quality score)\n"
"  -i  --input <path>\t\tLocation of input files [default: \"\"]\n"
"  -o  --output <path>\t\tLocation to create output files [default: \"\"]\n"
"  -e  --logfile <path>\t\tLocation of message output [default: none]\n"
"\t\t\t\t(Alternative to script redirect of error output)\n"
"  -l  --loglevel <level>\tLevel of message output [default: warning]\n"
"\t\t\t\t(none/fatal/error/warning/information/debug)\n"
"  -s  --simdata <header>\tOutput simulation data\n"
"  -w  --working\t\t\tOutput final working values\n"
"  -M  --M <filepath>\t\tPredetermined Crosstalk matrix file path\n"
"  -P  --P <filepath>\t\tPredetermined Phasing matrix file path\n"
"  -N  --N <filepath>\t\tPredetermined Noise matrix file path\n"
"  -S  --solver <solver>\t\tLinear equation solver to use for P matrix\n"
"\tOptions:\t\tls   - least squares, allow negatives [default]\n"
"\t\t\t\tzero - least squares then set negatives to zero\n"
"\t\t\t\tnnls - non-negative least squares\n"
"\n"
"  --help\t\t\tDisplay this help\n"
"  --licence\t\t\tDisplay AYB licence information\n"
"  --version\t\t\tDisplay AYB version information\n"
"\n"
