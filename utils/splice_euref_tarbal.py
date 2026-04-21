import glob
import os
import subprocess
from datetime import datetime
from dateutil import parser
from tqdm import tqdm
gfzrnx_path = os.path.split(os.getcwd())[0] + f'{os.sep}gfzrnx{os.sep}'

def main(date, indir, verbose=False):
    year = parser.parse(date).year
    mmdd = parser.parse(date).strftime('%m%d')
    folder = f'{indir}/{year}/{mmdd}/'
    files = sorted(glob.glob(folder + '/*.tar.gz'))[:-1]

    gfzrnx = glob.glob(gfzrnx_path + "gfzrnx*")[0]
    gkey = 'gzip'
    gflags = '-d -f -q'

    # UNZIP ALL
    if verbose:
        for i in tqdm(range(len(files)), desc="Unzipping"):
            try:
                print (f"Unzipping {files[i]}, {i+1}/{len(files)}")
                #Find the corresponding rx filename
                subprocess.call(f'{gkey} {gflags} {files[i]}', shell=True)
            except Exception as e:
                if verbose:
                    print (f"Error processing {files[i]}: {e}")
    else:
        for i in range(len(files)):
            try:
                # print (f"Unzipping {files[i]}, {i+1}/{len(files)}")
                #Find the corresponding rx filename
                subprocess.call(f'{gkey} {gflags} {files[i]}', shell=True)
            except Exception as e:
                pass

    files = sorted(glob.glob(folder + '/*.tar'))
    # TAR expand and splice all files
    for iff, f in enumerate(files):
        try:
            rx = os.path.basename(f).split('_')[0]
            print (f"Processing {f}, {iff+1}/{len(files)}")
            subdir = folder + f'/{rx}/'
            subprocess.call(f'mkdir -p {subdir}', shell=True)
            #TAR expand the tar file
            subprocess.call(f'tar -xf {f} -C {subdir}', shell=True)
            
            #Find all partial hatanaka compressed files
            subrx_file_path = sorted(glob.glob(f"{subdir}/{rx}*.crx"))
            #Decompress all hatanaka compressed files to rnx
            if verbose:
                for i in tqdm(range(len(subrx_file_path)), desc="Converting"):
                        if os.path.exists(subrx_file_path[i].replace('.crx', '.rnx')):
                            continue
                        subprocess.call(f'crx2rnx {subrx_file_path[i]}', shell=True)
            else:
                for i in range(len(subrx_file_path)):
                    if os.path.exists(subrx_file_path[i].replace('.crx', '.rnx')):
                        continue
                    subprocess.call(f'crx2rnx {subrx_file_path[i]}', shell=True)
            
            #Find all decompressed rnx files
            subrx_file_path = sorted(glob.glob(f"{subdir}/*.rnx"))
            finp = " ".join(subrx_file_path)
            #Prepare the command for GFZRNX to splice the rnx files and convert to crx
            fout = os.path.splitext(os.path.splitext(f)[0])[0] + '.rnx'
            command = f"{gfzrnx} -finp {finp} -fout {fout} -kv -q -satsys CEG -splice_memsave"
            t0 = datetime.now()
            subprocess.call(command, shell=True)
            if verbose:
                print (f"It took: {datetime.now() - t0} to complete processing {rx}")
            #Convert the spliced rnx file to hatanaka compressed crx to save 90% of the space
            subprocess.call(f'rnx2crx {fout} ', shell=True)

            #Remove all intermediate files
            subprocess.call(f'rm -rf {subdir}', shell=True)
            subprocess.call(f'rm -rf {fout}', shell=True)
            subprocess.call(f'rm -rf {f}', shell=True)
        except Exception as e:
            if verbose:
                print (f"Error processing {f}: {e}")

    return 0


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('date', help='2017-5-27, or 2017-251, or a comma-separed start,end dates', type=str)
    p.add_argument('indir', help='Input directory', type=str)
    p.add_argument('--verbose', help='Verbose output', action='store_true')
    
    args = p.parse_args()
    main(args.date, args.indir)
    