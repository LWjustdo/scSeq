# python3
# coding = utf-8
# Author: liang.wan, 2019-05-20
"""
Convert result.json to html report for DragonDropinn single cell sequence.

"""
import json
import os
import csv
#from functools import partial


class Report:
    """
    a class to generate a HTML report from a template html by use a result.json

    """
    def __init__(self,  prefix, template_html, result_json, out_dir=os.getcwd()):
        self.prefix = prefix
        self.templateHtml = template_html
        self.resultJson = result_json
        self.outDir = out_dir

    @staticmethod
    def tag(name, *content, cls=None, **attrs):
        """generate HTML tag"""
        if cls is not None:
            attrs["class"] = cls
        if attrs:
            attrs_str = "".join(" %s='%s'" % (attr, value)
                                for attr, value
                                in sorted(attrs.items()))
        else:
            attrs_str = ""
        if content:
            return "\n".join("<%s%s>%s</%s>" % (name, attrs_str, c, name)
                             for c in content)
        else:
            return "<%s%s />" % (name, attrs_str)


    @staticmethod
    def convert_table(infile):
        """convert a tsv or csv to a table tag string in html"""
        row = 0
        out = "<table>"
        postfix = os.path.splitext(infile)[-1]
        #print(postfix)
        if postfix == ".csv":
            with open(infile, "r") as f:
                f = csv.reader(f)
                for line in f:
                    row += 1
                    if row == 1:
                        out += "<tr>"
                        tmp = ["<td class='header_col'>" + i + "</td>" for i in
                               line]
                        out += "".join(tmp)
                        out += "</tr>"
                    elif 1 < row <= 10:
                        out += "<tr>"
                        tmp = ["<td>" + i + "</td>" for i in line]
                        out += "".join(tmp)
                        out += "</tr>"
                    else:
                        break
                out += "</table>"
        else:
            with open(infile, "r") as f:
                for line in f:
                    #if type(line) == str:
                    #    pass
                    #else:
                    #    line = line.decode('utf-8')
                    #print(line)
                    row += 1
                    line = line.strip().split("\t")
                    if row == 1:
                        out += "<tr>"
                        tmp = ["<td class='header_col'>" + i + "</td>" for i in
                               line]
                        out += "".join(tmp)
                        out += "</tr>"
                    elif 1 < row <= 10:
                        out += "<tr>"
                        tmp = ["<td>" + i + "</td>" for i in line]
                        out += "".join(tmp)
                        out += "</tr>"
                    else:
                        break
                out += "</table>"
        return out

    @staticmethod
    def choose_tag(akey, avalue):
        """choose which tag to write"""
        # if akey.endswith("图"):
        #    return picture(src=avalue, alt=akey, width=800)
        if akey.endswith("表"):
            return Report.convert_table(avalue)
        else:
            return avalue



    def read_template(self):
        """read template html return a string"""
        with open(self.templateHtml, "r") as f:
            template = f.read()
        if type(template) == str:
            return template
        else:
            return template.decode('utf-8')

    def write_html(self, html_string):
        """write output html to output directory"""
        outfile = os.path.join(self.outDir, self.prefix + "_scSeqReport.html")
        with open(outfile, "w") as f:
            f.write(html_string)
        return

    def read_and_check_json(self):
        """load result.json, the key of the json is same as the template html"""
        with open(self.resultJson, "r") as f:
            resultDic = json.load(f)
            projectInfoDic = resultDic["项目信息"]
            customerInfoDic = resultDic["客户信息"]
            sequenceInfoDic = resultDic["测序信息"]
            analysisInfoDic = resultDic["分析结果"]
            for k, v in projectInfoDic.items():
                projectInfoDic[k] = v if projectInfoDic[k] != "" else "NA"
            for k, v in customerInfoDic.items():
                customerInfoDic[k] = v if customerInfoDic[k] != "" else "NA"
            for k, v in sequenceInfoDic.items():
                sequenceInfoDic[k] = v if sequenceInfoDic[k] != "" else "NA"
            for k, v in analysisInfoDic.items():
                if os.path.exists(v):
                    pass
                else:
                    raise Exception("FileNotFoundError: " + k + " --> " + v +
                                    ", check and run again!")
        #print(projectInfoDic, customerInfoDic, sequenceInfoDic, analysisInfoDic)
        return projectInfoDic, customerInfoDic, sequenceInfoDic, analysisInfoDic

    def paste_html_mark(self):
        html = self.read_template()
        p, c, s, a = self.read_and_check_json()
        for k, v in p.items():
            html = html.replace("{%" + k + "%}", v)
        for k, v in c.items():
            html = html.replace("{%" + k + "%}", v)
        for k, v in s.items():
            html = html.replace("{%" + k + "%}", v)
        for k, v in a.items():
            tmp = Report.choose_tag(k, v)
            html = html.replace("{%" + k + "%}", tmp)
        return html











