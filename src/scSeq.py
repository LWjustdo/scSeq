"""
scSeq.py - Tools for DragonDropinn single cell sequence
=======================================================
@Author: liang.wan
@Date:2019-05-23

There are 9 functions:
 - summary4fastp
 - extract_head_seq
 - get_single_cell_read_count
 - map_tsv2_all
 - sequence_saturation
 - generate_report
 - grep_STAR_map_info
 - grep_seurat_cell_filter
 - generate_json
To get help on a specific func, type:
    scSeq.py <func> --help
To ues a specific func, type:
    scSeq.py <func> [func options]
"""
from scSeqReport import Report
from version import _version_
import scSeq_methods as sc
import sys

def main(argv=None):
    funcs = ["summary4fastp",
             "extract_head_seq",
             "get_single_cell_read_count",
             "map_tsv2_all",
             "sequence_saturation",
             "generate_report",
             "grep_STAR_map_info",
             "grep_seurat_cell_filter",
             "generate_json"]
    argv = sys.argv
    if len(argv) == 1 or argv[1] == "-h" or argv[1] == "--help":
        print(globals()["__doc__"])
        return
    elif len(argv) == 1 or argv[1] == "-v" or argv[1] == "--version":
        print("\nscSeq.py Version: %s.\n" % (_version_))
        return
    func = argv[1]
    if func in funcs:
        if func == "generate_report":
            if len(argv) <= 2 or argv[2] == "-h" or argv[2] == "--help":
                print("\nscSep.py generate_report <prefix> <template_html> "
                      "<result_json> <out_dir>\n")
                return
            else:
                try:
                    report = Report(argv[2], argv[3], argv[4], argv[5])
                    html = report.paste_html_mark()
                    report.write_html(html)
                except Exception as e:
                    print(e)
                    print("Arguments error: 'scSeq.py <func> -h for help'")
                    return
        else:
            function = getattr(sc, func)
            if len(argv) <= 2 or argv[2] == "-h" or argv[2] == "--help":
                print(function.__doc__)
                return
            else:
                del argv[0]
                del argv[0]
                function(*argv)
                #try:
                #    function(*argv)
                #except:
                #    print("\nArguments error: 'scSeq.py <func> -h for help'\n")
                #    return
    else:
        print("\n'%s' is not defined, type 'scSeq.py -h for help'!\n" %(func))
        return


if __name__ == "__main__":
    sys.exit(main())

