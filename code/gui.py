
import tkinter as tk
from tkinter.filedialog import askopenfilename
import subprocess, os

f1 = None # RNA/DNA file path
f2 = None # ligands file path

def browsefunc(label, index):
    """ Browse for input files """
    filepath = askopenfilename()
    label.config(text = filepath.split("/")[-1])

    global f1, f2
    if index == '1':
        f1 = filepath
    else:
        f2 = filepath

def run_script():
    """ Run fingeRNAt.py after pressing Submit button """

    global f1, f2
    if f1 is None:
        lbl5.config(text = "Choose input RNA/DNA file first!", fg = "red")

    else:
        dha = "" if option2.get() == "No" else "-dha"
        wrapper = "" if option3.get() == "None" else "-wrapper %s" %option3.get()
        addH = "" if option4.get() == "No" else "-addH %s" %option4.get()
        h2o = "" if optionwater.get() == "No" else "-h2o"
        detail = "" if optiondetail.get() == "No" else "-detail"
        yaml = "" if len(yaml_entry.get()) == 0 else "-custom %s" %yaml_entry.get()
        out = "" if len(entry.get()) == 0 else "-o %s" %entry.get()

        if f2 is None:
            subprocess.call('python %s/fingeRNAt.py -r %s -f %s %s %s %s %s %s %s' % (os.path.dirname(os.path.realpath(__file__)), f1, option.get(), dha, wrapper, addH, out, h2o, detail), shell = True)
        else:
            subprocess.call('python %s/fingeRNAt.py -r %s -l %s -f %s %s %s %s %s %s %s %s' % (os.path.dirname(os.path.realpath(__file__)), f1, f2, option.get(), dha, wrapper, addH, out, h2o, detail, yaml), shell = True)

root = tk.Tk()
root.title("fingeRNAt")
root.geometry("575x900")
root.resizable(width = False, height = True)

# Intro

app_intro = tk.Frame(root,bd=16, relief='sunken')
app_intro.grid(row = 0, column = 0)

intro =  tk.Label(app_intro, text = "Welcome to fingeRNAt program!")
intro.grid(row = 0, column = 0, sticky = tk.W + tk.E, pady = 5, columnspan = 2)

intro2 = tk.Label(app_intro, text = "Choose your inputs and parametres to calculate Structural Interaction Fingerprint (SIFt)")
intro2.grid(row = 1, column = 0, sticky = tk.W + tk.E, columnspan = 2)

# Input files choice

app1 = tk.Frame(root)
app1.grid(row = 1, column=  0, pady = 7)

lbl1 = tk.Label(app1, text = "Choose path to RNA/DNA file")
lbl1.grid(row = 0, column = 0, sticky = tk.W )

input1_choice = tk.Label(app1, text="No file chosen", fg = "green")
input1_choice.grid(row = 2, column = 0)

browse_bttn1 = tk.Button(app1, text = "Browse...", highlightbackground="grey", fg = "white", command = lambda: browsefunc(input1_choice,'1'))
browse_bttn1.grid(row = 1, column = 0)

lbl2 = tk.Label(app1, text = "Choose path to ligands file")
lbl2.grid(row = 0, column = 1, sticky = tk.W )

input2_choice = tk.Label(app1, text="No file chosen", fg = "green")
input2_choice.grid(row = 2 ,column = 1)

browse_bttn2 = tk.Button(app1, text = "Browse...", highlightbackground = "grey", fg = "white", command = lambda: browsefunc(input2_choice,'2'))
browse_bttn2.grid(row = 1, column = 1)

# Choose SIFt type

app2 = tk.Frame(root)
app2.grid(row = 2, column = 0, pady = 7)

lbl3 = tk.Label(app2, text = "Choose Structural Interaction Fingerprint (SIFt) type")
lbl3.grid(row = 0, column = 0, columnspan = 2, sticky = tk.W + tk.E)
option = tk.StringVar(None,"FULL")

choicebttn1 = tk.Radiobutton(app2, variable = option, text = "FULL", value = "FULL")
choicebttn1.config(indicatoron = 0, bd = 4, width = 10)
choicebttn1.grid(row = 1, column = 0, columnspan=20)

choicebttn2 = tk.Radiobutton(app2, variable = option, text = "SIMPLE", value = "SIMPLE")
choicebttn2.config(indicatoron = 0, bd = 4, width = 10)
choicebttn2.grid(row = 2, column = 0, sticky = tk.E)

choicebttn3 = tk.Radiobutton(app2, variable = option, text = "PBS", value = "PBS")
choicebttn3.config(indicatoron = 0, bd = 4, width = 10)
choicebttn3.grid(row=2, column=1, sticky = tk.W)


# How to add hydrogens to ligands structures

app4 = tk.Frame(root)
app4.grid(row = 3, column = 0, pady = 7)

lbl4 = tk.Label(app4, text = "Choose module to add hydrogens to ligands")
lbl4.grid(row = 0, column = 0, columnspan = 2, sticky = tk.W + tk.E)

option4 = tk.StringVar(None,"OpenBabel")

choicebttn4_1 = tk.Radiobutton(app4, variable = option4, text="OpenBabel", value = 'OpenBabel')
choicebttn4_1.config(indicatoron = 0, bd = 4, width =10)
choicebttn4_1.grid(row = 1, column = 0, columnspan=20)

choicebttn4_2 = tk.Radiobutton(app4, variable = option4, text="RDKit", value = 'RDKit')
choicebttn4_2.config(indicatoron = 0, bd = 4, width = 10)
choicebttn4_2.grid(row = 2, column = 0, sticky = tk.E)

choicebttn4_2 = tk.Radiobutton(app4, variable = option4, text="None", value = 'None')
choicebttn4_2.config(indicatoron = 0, bd = 4, width = 10)
choicebttn4_2.grid(row = 2, column = 1, sticky = tk.W)


# Consider dha?
app3_1 = tk.Frame(root)
app3_1.grid(row = 4, column = 0, pady = 7)

lbl3_1 = tk.Label(app3_1, text = "Consider Donor - Hydrogen - Acceptor angle when detecting hydrogen bonds?")
lbl3_1.grid(row = 0, column = 0, columnspan = 2, sticky = tk.W + tk.E)
lbl3_1 = tk.Label(app3_1, text = "(Applies only to SIFt type FULL)")
lbl3_1.grid(row = 1, column = 0, columnspan = 2, sticky = tk.W + tk.E)
lbl3_1.config(font=(None, 12))

option2 = tk.StringVar(None,"No")

choicebttn3_11 = tk.Radiobutton(app3_1, variable = option2, text="No", value = 'No')
choicebttn3_11.config(indicatoron = 0, bd = 4, width =10)
choicebttn3_11.grid(row = 2, column = 0, sticky = tk.E)

choicebttn3_12 = tk.Radiobutton(app3_1, variable = option2, text="Yes", value = 'Yes')
choicebttn3_12.config(indicatoron = 0, bd = 4, width = 10)
choicebttn3_12.grid(row = 2, column = 1, sticky = tk.W)


# Consider H2O?

app3_water = tk.Frame(root)
app3_water.grid(row = 5, column = 0, pady = 7)

lbl3_water = tk.Label(app3_water, text = "Consider water-mediated interactions?")
lbl3_water.grid(row = 0, column = 0, columnspan = 2, sticky = tk.W + tk.E)
lbl3_water = tk.Label(app3_water, text = "(Applies only to SIFt type FULL)")
lbl3_water.grid(row = 1, column = 0, columnspan = 2, sticky = tk.W + tk.E)
lbl3_water.config(font=(None, 12))

optionwater = tk.StringVar(None,"No")

choicebttn_water_1 = tk.Radiobutton(app3_water, variable = optionwater, text="No", value = 'No')
choicebttn_water_1.config(indicatoron = 0, bd = 4, width =10)
choicebttn_water_1.grid(row = 2, column = 0, sticky = tk.E)

choicebttn_water_12 = tk.Radiobutton(app3_water, variable = optionwater, text="Yes", value = 'Yes')
choicebttn_water_12.config(indicatoron = 0, bd = 4, width = 10)
choicebttn_water_12.grid(row = 2, column = 1, sticky = tk.W)

# Wrapper

app3_2 = tk.Frame(root)
app3_2.grid(row = 6, column = 0, pady = 7)

lbl3_2 = tk.Label(app3_2, text = "Choose wrapper")
lbl3_2.grid(row = 0, column = 0, columnspan = 2, sticky = tk.W + tk.E)

option3 = tk.StringVar(None,"None")

choicebttn3_21 = tk.Radiobutton(app3_2, variable = option3, text="None", value = 'None')
choicebttn3_21.config(indicatoron = 0, bd = 4, width =10)
choicebttn3_21.grid(row = 1, column = 0, sticky = tk.E)

choicebttn3_22 = tk.Radiobutton(app3_2, variable = option3, text="ACUG", value = 'ACUG')
choicebttn3_22.config(indicatoron = 0, bd = 4, width = 10)
choicebttn3_22.grid(row = 1, column = 1, sticky = tk.W)

choicebttn3_23 = tk.Radiobutton(app3_2, variable = option3, text="PuPy", value = 'PuPy')
choicebttn3_23.config(indicatoron = 0, bd = 4, width = 10)
choicebttn3_23.grid(row = 2, column = 0, sticky = tk.W)

choicebttn3_23 = tk.Radiobutton(app3_2, variable = option3, text="Counter", value = 'Counter')
choicebttn3_23.config(indicatoron = 0, bd = 4, width = 10)
choicebttn3_23.grid(row = 2, column = 1, sticky = tk.W)


# Save -detail?

app_detail = tk.Frame(root)
app_detail.grid(row = 7, column = 0, pady = 7)

lbl3_detail = tk.Label(app_detail, text = "Save details about interactions?")
lbl3_detail.grid(row = 0, column = 0, columnspan = 2, sticky = tk.W + tk.E)

optiondetail = tk.StringVar(None,"No")

choicebttn_detail_1 = tk.Radiobutton(app_detail, variable = optiondetail, text="No", value = 'No')
choicebttn_detail_1.config(indicatoron = 0, bd = 4, width =10)
choicebttn_detail_1.grid(row = 2, column = 0, sticky = tk.E)

choicebttn_detail_12 = tk.Radiobutton(app_detail, variable = optiondetail, text="Yes", value = 'Yes')
choicebttn_detail_12.config(indicatoron = 0, bd = 4, width = 10)
choicebttn_detail_12.grid(row = 2, column = 1, sticky = tk.W)

# Optional filename

app3 = tk.Frame(root)
app3.grid(row = 8, column = 0, pady = 7)

lbl4 = tk.Label(app3, text = "Optional output file path")
lbl4.grid(row = 0, column = 0, columnspan = 2, sticky = tk.W + tk.E)
entry = tk.Entry(app3,width=50)
entry.grid(row = 1, column = 0, columnspan = 2, sticky = tk.W + tk.E)

# Optional run fingerDISt directly on generated SIFts outputs

app_yaml = tk.Frame(root)
app_yaml.grid(row = 9, column = 0, pady = 5)

lbl_yaml = tk.Label(app_yaml, text = "Optional path to YAML file with user-defined interactions")
lbl_yaml.grid(row = 0, column = 0, columnspan = 2, sticky = tk.W + tk.E)
lbl_yaml = tk.Label(app_yaml, text = "(Applies only to SIFt type FULL)")
lbl_yaml.grid(row = 1, column = 0, columnspan = 2, sticky = tk.W + tk.E)
lbl_yaml.config(font=(None, 12))

yaml_entry = tk.Entry(app_yaml,width=50)
yaml_entry.grid(row = 2, column = 0, columnspan = 2, sticky = tk.W + tk.E)


# Submit button

app5 = tk.Frame(root)
app5.grid(row = 10, column = 0)

lbl5 = tk.Label(app5, text = "")
lbl5.grid(row = 0, column = 0)

bttn5 = tk.Button(app5, text = "Submit", state = "normal", width = 10, height = 2, highlightbackground = "blue", fg = "white", command = run_script)
bttn5.grid(row = 1, column = 0)

root.mainloop()
