import sys
from matplotlib_venn import venn2
from matplotlib import pyplot

output_fn = sys.argv[1]
output_image_fn = sys.argv[2]
input_list1_fn = sys.argv[3]
input_list2_fn = sys.argv[4]

with open(input_list1_fn) as f:
    offtarget_list = f.read().splitlines()
with open(input_list2_fn) as f:
    dTPM_list = f.read().splitlines()
offtarget_list_set = set(offtarget_list)
dTPM_list_set = set(dTPM_list)

offtarget_vs_dTPM_list_set_intersection = offtarget_list_set & dTPM_list_set # 積集合
uni_offtarget_difference = offtarget_list_set - offtarget_vs_dTPM_list_set_intersection # 差集合（offtarget）
uni_dTPM_difference = dTPM_list_set - offtarget_vs_dTPM_list_set_intersection # 差集合（tpm）

with open(output_fn, mode='w', encoding='utf-8') as fp:
    for line in list(offtarget_vs_dTPM_list_set_intersection):
        # write each item on a new line
        fp.write("%s\n" % line)

venn2(subsets=(len(uni_offtarget_difference),len(uni_dTPM_difference),len(offtarget_vs_dTPM_list_set_intersection)))

pyplot.savefig(output_image_fn, format="tiff", dpi=350)
