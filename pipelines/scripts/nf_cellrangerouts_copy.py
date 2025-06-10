import os
import shutil
import argparse

parser = argparse.ArgumentParser(description='Process files based on parameters')
parser.add_argument('--datlabel', help='Dataset name', required=True)
parser.add_argument('--cellrangermulti', help='Whether Cell Ranger output is multi (0 or 1)', type=int, choices=[0, 1], required=True)
parser.add_argument('--cellrangeroutpath', help='Path to Cell Ranger output directory', required=True)
args = parser.parse_args()

# Handle dataset directory creation
if not os.path.isdir(args.datlabel):
    os.makedirs(args.datlabel)

cellrangermulti = bool(args.cellrangermulti)


# Handle multi condition
if args.cellrangermulti:
    for sample in os.listdir(os.path.join(args.cellrangeroutpath, "outs/per_sample_outs/")):
        sample_path = os.path.join(args.cellrangeroutpath, "outs/per_sample_outs", sample)
        count_path = os.path.join(sample_path, "count/sample_filtered_feature_bc_matrix/")
        if os.path.isdir(count_path):
            for file in os.listdir(count_path):
                src = os.path.join(count_path, file)
                dest = os.path.join(args.datlabel, f"{args.datlabel}_{sample}_{file}")
                shutil.copy(src, dest)
else:
    input_files = os.listdir(os.path.join(args.cellrangeroutpath, "outs/filtered_feature_bc_matrix/"))
    for file in input_files:
        src = os.path.join(args.cellrangeroutpath, "outs/filtered_feature_bc_matrix", file)
        dest = os.path.join(args.datlabel, f"{args.datlabel}_{file}")
        shutil.copy(src, dest)
