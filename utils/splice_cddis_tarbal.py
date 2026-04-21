import glob
import os
import subprocess
from datetime import datetime
from dateutil import parser
from tqdm import tqdm
gfzrnx_path = os.path.split(os.getcwd())[0] + f'{os.sep}gfzrnx{os.sep}'

def main(date, indir, verbose=False):
    year = parser.parse(date).year
    doy = parser.parse(date).timetuple().tm_yday
    YY = str(year)[-2:] 
    mmdd = parser.parse(date).strftime('%m%d')
    folder = f'{indir}/{year}/{mmdd}/'
    files = sorted(glob.glob(folder + '/*.tar'))
    sub_path = f'/gnss/data/highrate/{year}/{doy}/{YY}d/'

    gfzrnx = glob.glob(gfzrnx_path + "gfzrnx*")[0]
    gkey = 'gzip'
    gflags = '-d -f -q'

    for iff, f in enumerate(files):
        try:
            print (f"Processing {f}, {iff+1}/{len(files)}")
            #Find the corresponding rx filename
            rx = os.path.basename(f).split('_')[0]
            #TAR expand the tar file
            subprocess.call(f'tar -xf {f} -C {folder}', shell=True)
            #Find all partial files in subdirectories
            subrx_file_path_gz = sorted(glob.glob(f"{folder + sub_path}/*/{rx}*.crx.gz"))
            #Unzip all partial files
            if verbose:
                for i in tqdm(range(len(subrx_file_path_gz)), desc="Unzipping"):
                    subprocess.call(f'{gkey} {gflags} {subrx_file_path_gz[i]}', shell=True)
            else:
                for i in range(len(subrx_file_path_gz)):
                    subprocess.call(f'{gkey} {gflags} {subrx_file_path_gz[i]}', shell=True)
            #Find all new unzipped hatanaka compressed files
            subrx_file_path = sorted(glob.glob(f"{folder + sub_path}/*/{rx}*.crx"))
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
            #Remove all hatanaka partial files
            subprocess.call(f'rm -rf {" ".join(subrx_file_path)}', shell=True)
            #Find all decompressed rnx files
            subrx_file_path = sorted(glob.glob(f"{folder + sub_path}/*/{rx}*.rnx"))
            finp = " ".join(subrx_file_path)
            #Prepare the command for GFZRNX to splice the rnx files and convert to crx
            fout = os.path.splitext(os.path.splitext(f)[0])[0] + '.rnx'
            command = f"{gfzrnx} -finp {finp} -fout {fout} -kv -q -satsys CEG -splice_memsave"
            t0 = datetime.now()
            subprocess.call(command, shell=True)
            print (f"It took: {datetime.now() - t0} to complete processing {rx}")
            #Convert the spliced rnx file to hatanaka compressed crx to save 90% of the space
            subprocess.call(f'rnx2crx {fout} ', shell=True)
            #Remove all intermediate files
            subprocess.call(f'rm -rf {finp}*', shell=True)
            subprocess.call(f'rm -rf {fout}', shell=True)
            subprocess.call(f'rm -rf {f}', shell=True)
        except Exception as e:
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
    