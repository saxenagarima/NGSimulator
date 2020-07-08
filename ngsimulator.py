from Tkinter import *
import tkFileDialog as filedialog
import os
from PIL import ImageTk, Image
from os.path import basename
from tkinter.font import Font
class MyFrame:
	def __init__(self, root):
		self.master=root
		self.master.title("In Silico Sequencer")

		##UHTL logo
		'''
		load=Image.open("UHTL.png")
		render=ImageTk.PhotoImage(load)
		img=Label(root, image=render)
		img.image=render
		img.place(x=0,y=0,height=50, width=200)
		'''

		version_info=Text(root)
		version_info.place(x=400, y=5,height=50, width=210)
		version_info.insert(END, "In Silico Sequencer (v1.0.0)"+"\n"+"\n"+"Developed by:Mehul Jani")
		version_info.configure(state='disabled', bg=root["bg"])

		aa=10 
		bb=65
		##Upload Reference file
		buttonr = Button(self.master, text="Browse", command=self.load_file)
		buttonr.place(x=aa,y=bb+20,height=25, width=75)	
	
		self.filename=Entry(self.master)
		self.filename.place(x=aa+75,y=bb+20,height=25, width=525)
		
		reflab=Label(self.master, text="Please upload a reference file in fasta format")
		reflab.place(x=aa, y=bb)
		
		##Upload Bedfile
		buttonb = Button(self.master, text="Browse", command=self.load_filebed)
		buttonb.place(x=aa,y=bb+80,height=25, width=75)	
	
		self.filename1=Entry(self.master)
		self.filename1.place(x=85,y=bb+80,height=25, width=525)
		
		bedlab=Label(self.master, text="Please upload a bedfile in tsv format")
		bedlab.place(x=aa, y=bb+60)

		##Entryboxes
		label2=Label(self.master, text="Read length")
		label2.place(x=aa, y=bb+120)

		self.Len_entry=Entry(self.master)
		self.Len_entry.place(x=aa+150, y=bb+115, height=25, width=50)
		self.Len_entry.insert(END, '100')
		
		
		coverage=Label(self.master, text="Coverage")
		coverage.place(x=aa, y=bb+150)

		self.cov_entry=Entry(self.master)
		self.cov_entry.place(x=aa+150, y=bb+145, height=25, width=50)
		self.cov_entry.insert(END, '100')		
		

		label3=Label(self.master, text="SNP error rate")
		label3.place(x=aa, y=245 )

		self.ser_entry=Entry(self.master)
		self.ser_entry.place(x=aa+150, y=bb+180, height=25, width=50)
		self.ser_entry.insert(END, '0.5')
		

		label4=Label(self.master, text="Insertion error rate")
		label4.place(x=aa+290, y=bb+120 )

		self.ier_entry=Entry(self.master)
		self.ier_entry.place(x=aa+450, y=bb+115, height=25, width=50)
		self.ier_entry.insert(END, '0.5')
		
		
		label5=Label(self.master, text="Deletion error rate")
		label5.place(x=aa+290, y=bb+150 )

		self.der_entry=Entry(self.master)
		self.der_entry.place(x=aa+450, y=bb+145, height=25, width=50)
		self.der_entry.insert(END, '0.5')	
		

		label6=Label(self.master, text="Region(s)")
		label6.place(x=aa+290, y=bb+180 )

		self.reg_entry=Entry(self.master)
		self.reg_entry.place(x=aa+450, y=bb+175, height=25, width=150)
		self.reg_entry.insert(END, '200:500 800:1000')

		##Run button
		runbutton=Button(self.master, text="Run", bg="pale green", command=self.gen_reads, width=10)
		runbutton.place(x=aa,y=bb+225,height=50, width=75)

	def load_file(self):
		self.file_selected = filedialog.askopenfilename(filetypes =(("Fasta", "*.fasta"),("All Files","*.*")))
		self.filename.delete(0,END)
		self.filename.insert(0,self.file_selected)
		if self.file_selected:
			print self.file_selected, "file selected"
		else:
			print "Please select a fastq file"

	def load_filebed(self):
		self.file1_selected = filedialog.askopenfilename(filetypes =(("Tab separated", "*.tsv"),("All Files","*.*")))
		self.filename1.delete(0,END)
		self.filename1.insert(0,self.file1_selected)
		if self.file1_selected:
			print self.file1_selected, "file selected"
		else:
			print "Please select a tsv file"

	def gen_reads(self):
		print self.Len_entry.get(),self.cov_entry.get(),self.ser_entry.get(),self.der_entry.get(),self.ier_entry.get(),self.reg_entry.get(),self.file_selected,self.file1_selected
		#print os.path.splitext(self.file_selected)[0]+".fastq",os.path.splitext(self.file1_selected)[0]+".out"
		outname=os.path.splitext(os.path.basename(self.file_selected))[0]
		os.system("python readgen.py %s %s %s -len %s -cov %s -ser %s -der %s -ier %s -hsr %s > %s" %(self.file_selected, self.file1_selected,os.path.splitext(self.file_selected)[0]+".fastq",self.Len_entry.get(), self.cov_entry.get(),self.ser_entry.get(),self.der_entry.get(),self.ier_entry.get(), self.reg_entry.get(), os.path.splitext(self.file_selected)[0]+".out"))
		self.label6=Label(self.master, text="Please check %s and %s files!!"%(outname+".fastq",outname+".out"))
		self.label6.place(x=110, y=305)
		self.flash()
		
	def flash(self):
		try:
			self.prv
		except:
			fg = self.label6.cget("foreground")
		else:
			fg=self.prv
			if fg is "black":
				fg="gray64"
			else:
				fg="black"
		
		self.label6.configure(foreground=fg)
		self.prv=fg
		self.label6.after(1000, self.flash)

def main():
	root=Tk()
	MyFrame(root)
	root.geometry("650x360")
	root.resizable(0,0)
	root.mainloop()


if __name__ == '__main__':
    main()  
