HOW TO RUN THE PROGRAM

install Python: https://www.python.org/downloads/
run "setup.py" with python (you only need to do these steps once on every computer you want to run the program on)

now you can run the program "DIANN_to_finished_output.py" with python

WHAT THE PROGRAM DOES

The program takes DIA-NN output (proteingroups), as a comma seperated txt file
and returns filtered, log2 transformed values, and if all the values in a group are absent, they are replaced with 1

The program will promt you for 2 file names
The first input needs to be the DIANN output, as a csv, with ";" seperator
The second input needs to be an (empty) .txt file, into which the new will be added

You can now import all your data into perseus and start visualizing ;)

GROUP NOMENCLATURE 

The program groups the columns based on the regular expression "(d\d+)_\d+", which has to appear somewhere in the column name
it looks in the column name for the phrase dx_x, with x being any amount of decimals
and groups the columns based on the first decimal

Example:

	T: Protein Groups | d0_1 | d0_2 | d0_3 | d1_1 | d1_2 | d1_3 |

	->

	T: Protein Groups | d0_1 | d0_2 | d0_3 | d1_1 | d1_2 | d1_3 |
	groups            | d0   | d0   | d0   | d1   | d1   | d1   |


The expression just has to be somewhere in the column title, so file names also work
If using file names, make sure the regular expression doesnt appear anywhere else in the file name

Example 2: 

	T: Protein Groups | C://.../d0_1 | C://.../d0_2 | C://.../d0_3 | C://.../d1_1 | C://.../d1_2 | C://.../d1_3 |

	->

	T:  Protein Groups | C://.../d0_1 | C://.../d0_2 | C://.../d0_3 | C://.../d1_1 | C://.../d1_2 | C://.../d1_3 |
	groups             | d0           | d0           | d0           | d1           | d1           | d1           |


CHANGING SETTINGS

open "DIANN_to_finished_output.py" with a text editor:

to change valid percentage: 

line 10: filt = get_invalid(df1, percent=70) 
	-> change 70 to the desired percentage

to change impute settings:

line 17: df2 = impute(df1, width=0.3, downshift=1.8) 
	-> change width, downshift to desired values


