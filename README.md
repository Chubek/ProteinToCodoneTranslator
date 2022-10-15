# Pal2Nal

# Note
Run `python3 pal2nal.py -h` to see a detailed guide. PLEASE do it before you attempt to make a run. However, with the test files provided, you can run the script without any arguments to quickly see and time the script.

# How to Run
```
$ sudo apt-get install python3-venv #in case you don't have it
$ python3 -m venv env
$ source env/bin/activate
$ pip3 install -r requirements.txt
$ python3 pal2nal.py -h
# possible commands
$ python3 pal2nal.py -t genetic_table.json
$ python3 pal2nal.py -i 2 #translation table ID
$ python3 pal2nal.py -t https://file.share/my_genetic_table.json -aaf https://file.share/my_fasta.aa.fa -ntf https://file.share/my_fasta.nt.fa -tr my_fasta_result.fa #accepts links alternatively
$ python3 pal2nal.py -t (123IsuiiksK, genetic_table_google_drive.json) #accepts Google Drive links in this format. Google Drive ID is in its URL.
$ python3 pal2nal.py -i 4 -ntf https://file.share/my_nt_1.nt.fa my_nt_2.nt.py -aaf my_aa_1.aa.fa my_aa_2.aa.fa -tr my_result_1.nt.fa my_result_2.nt.fa #accepts multiple args
$ python3 pal2nal.py -i 4 -ntf my_nt_1.nt.fa my_nt_2.nt.py -aaf my_aa_1.aa.fa my_aa_2.aa.fa #if you don't pass a target it will save files iteratively with their index attached e.g. my_nt_2_result_2.nt.fa
$ python3 pal2nal.py -i 4 -ntf my_nt_1.nt.fa my_nt_2.nt.py -aaf my_aa_1.aa.fa my_aa_2.aa.fa -tr accumulative_result.nt.fa #if you pass one target with multiple inputs it will accumulate them into one file
```

Note: the multiple file functionality is untested as of beta-1.

Contact: Chubak#7400