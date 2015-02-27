#!/usr/bin/python3.4

######################################################
#Program Name:Robert Battocletti
#Date Published: 26FEB2014
#Designed to evaluate fastq files for microsatellites
######################################################
import re
from itertools import product
import configparser
from glob import glob
from sys import exit
from os import system


class ConfRead:
    def __init__(self):
        self.argvr = []
        self.config_issues = []
        config = configparser.ConfigParser()
        config.readfp(open(r'microsatellites.conf'))
        self.return_tetra = str(config.get('user_defined', 'return_tetra')).lower()
        self.min_adjacent_tetra = int(config.get('user_defined', 'min_adjacent_tetra'))
        self.min_adjacent_tri = int(config.get('user_defined', 'min_adjacent_tri'))
        self.min_adjacent_di = int(config.get('user_defined', 'min_adjacent_di'))
        self.lead_filter = int(config.get('user_defined', 'lead_filter'))
        self.rear_filter = int(config.get('user_defined', 'rear_filter'))
        self.project_name = config.get('user_defined', 'project_name')
        self.permutation_elements = (2, 3, 4)
        self.filter_tetra_min_adjacent = int(config.get('user_defined', 'filter_tetra_min_adjacent'))
        self.filter_tetra_min_nonadjacent = int(config.get('user_defined', 'filter_tetra_min_non-adjacent'))
        self.min_non_adjacent_tetra = int(config.get('user_defined', 'min_non_adjacent_tetra'))
        self.filter_tri_min_adjacent = int(config.get('user_defined', 'filter_tri_min_adjacent'))
        self.filter_tri_min_nonadjacent = int(config.get('user_defined', 'filter_tri_min_non-adjacent'))
        self.filter_di_min_adjacent = int(config.get('user_defined', 'filter_di_min_adjacent'))
        self.filter_di_min_nonadjacent = int(config.get('user_defined', 'filter_di_min_non-adjacent'))
        self.return_di = str(config.get('user_defined', 'return_di')).lower()
        self.return_tri = str(config.get('user_defined', 'return_tri')).lower()
        self.return_potential = str(config.get('user_defined', 'return_potential')).lower()
        self.process_directory = int(config.get('computing', 'process_directory'))
        self.print_result_to_screen = int(config.get('computing', 'print_result_to_screen'))
        self.max_non_adjacent_filter_tetra = int(config.get('user_defined', 'max_non_adjacent_filter_tetra'))
        self.read_average_length = int(151)
        self.write_database = "no"
        self.result_file = self.project_name + "_microsatellites.txt"
        self.database = self.project_name + ".db"
        self.result_csv = self.project_name + "_results.csv"
        self.potential_results = self.project_name + "_potential.txt"

        #self.valid() need to fix them
        self.read_directory()
        self.output_files()

    def __repr__(self):
        return repr((self.project_name, self.fn))

    def current_conf(self):
        cconf = '\n----------------\nCurrent configuration summary:\n----------------\nProject Name: {0}\nProcess folder:{7}\n--------------\nReturn tetra: {1}\nMin_adjacent tetra: {2}\nReturn tri: {3}\nMin_adjacent tri: {4}\nReturn di: {5}\nMin adjacent di: {6}\n--------------\nLead filter area: {8}\nRear filter area: {9}\nFilter tetra min non-adjacent: {10}\nFilter tetra min adjacent: {11}\nFilter tri min non-adjacent: {12}\nfilter tri min adjacent: {13}\nfilter di min non-adjacent:{14}\nfilter di min adjacent: {15}\n----------------\nPossible long repeats\n----------------\nReturn postential large repeats: {16}\nmin_non_adjacent_tetra: {17}\nmax_non_adjacent_filter_tetra: {18}'
        return cconf.format(self.project_name, self.return_tetra, self.min_adjacent_tetra, self.return_tri, self.min_adjacent_tri, self.return_di, self.min_adjacent_di, self.process_directory, self.lead_filter,
                            self.rear_filter, self.filter_tetra_min_nonadjacent, self.filter_tetra_min_adjacent, self.filter_tri_min_nonadjacent, self.filter_tri_min_adjacent, self.filter_di_min_nonadjacent, self.filter_di_min_adjacent, self.return_potential, self.min_non_adjacent_tetra, self.max_non_adjacent_filter_tetra)

    def output_files(self, num = 0):
        try:
            with open(self.result_file, "a") as result_file:
                result_file.write(self.current_conf())
        except Exception as exception:
            print(exception)
        if num == 0:
            if self.write_database == "yes":
                try:
                    with open(self.database, "a"):
                        pass
                except Exception as exception:
                    print(exception)
            try:
                with open(self.result_csv, "a") as csv_file:
                    csv_file.write("Name,Sequence\n")
            except Exception as exception:
                print(exception)
            if self.return_potential == "yes":
                try:
                    with open(self.potential_results, "a"):
                        pass
                except Exception as exception:
                    print(exception)

    def read_directory(self):
        fn = []
        if self.process_directory == 0:
            fn.append(input("Enter the file you would like to process :> "))
        elif self.process_directory == 1:
            switch = 0
            while not switch:
                fn_tmp = []
                files = glob("*fastq")
                for file in files:
                    print(file)
                    fn_tmp.append(file)
                print("\nI will process the following list is this correct enter(Y or press enter)? enter (N) to exit")
                prompter = input(str(":> ")).lower()
                if prompter == "n":
                    input("\nUpdate your working folder with the files you want to process and then press enter.\n")
                else:
                    fn = fn_tmp; switch = 1
        else:
            print("I believe your dna.conf is missing a vale for process_directory")
        self.fn = fn

    def valid(self):
        # Need to fix the first 3 validation statements
        filter_sum = self.lead_filter + self.rear_filter
        self.min_adjacent_tetra = 0 if self.return_tetra == 'yes' and self.min_adjacent_tetra < 151 - filter_sum else self.config_issues.append(
            """-return_tetra or min_adjacent_tetra not properly set.\n\t1) Ensure return_tetra is set to yes or no.\t\n2) Ensure min_adjacent_tetra + lead_filter + rear_filter is less than 151""")
        self.min_adjacent_tri = 0 if self.return_tri == 'yes' and self.min_adjacent_tri < 151 - filter_sum else self.config_issues.append(
            """-return_tri or min_adjacent_tri not properly set.\n\t1) Ensure return_tri is set to yes or no.\t\n2) Ensure min_adjacent_tri + lead_filter + rear_filter is less than 151""")
        self.min_adjacent_di = 0 if self.return_di == 'no' and self.min_adjacent_di < 151 - filter_sum else self.config_issues.append(
            """-return_di or min_adjacent_di not properly set.\n\t1) Ensure return_di is set to yes or no.\t\n2) Ensure min_adjacent_di + lead_filter + rear_filter is less than 151""")
        if filter_sum > int(120): self.config_issues.append(
            "The filter region is set properly: lead {}, rear {}.\n\t1) Recommend a value between 0-50.".format(
                self.lead_filter, self.rear_filter))
        if self.return_potential != 'no': self.config_issues.append("-return_potential is not properly set.")
        if self.return_potential == "yes" and self.min_non_adjacent_tetra == 0: self.config_issues.append(
            "- return_potential or min_non_adjacent_tetra not set properly")
        #self.config_issues.append("-process_directory is not properly set.") if self.process_directory != 0 or if self.process_directory != 1
        if self.config_issues:
            for issue in self.config_issues:
                print(issue)
                print("\nPlease fix the issues above and restart the program.")
                exit()

class Matcher:
    def __init__(self, argr):
        self.conf = argr
        self.permutation_generator()
        self.matrix()
        self.current_id = ''
        self.current_fn = ''
        self.current_header = ''
        self.current_read = ''
        self.current_read_leadfilter = self.current_read[:self.conf.lead_filter]
        self.current_read_rearfilter = self.current_read[-self.conf.rear_filter:]
        self.current_results = {}
        self.current_pos_results = {}
        self.current_counter = int(0)
        self.current_result = 0
        self.current_dirty = 0
        self.current_match_occur_in_read = 0
        self.current_match_occur_outside_filter = 0

    def cleaner(self):
        self.current_fn = ''
        self.current_header = ''
        self.current_read = ''
        self.current_id = ''
        self.current_results = {}
        self.current_pos_results = {}
        self.current_result = 0
        self.current_dirty = 0

    def intro_clear(self):
        lambda : os.system('clear')

    def intro(self, num):
        if num == 0:
            print("\n\n")
            input("Ensure that the dna.py and dna.conf are in the same folder. \npress enter to continue\n\n")
            self.intro_clear()
        elif num == 1:
            print("This process could take up to an one hour per 16GB file.")
            print("When complete, there will three different output files:")
            print("one database, one text file with results and a csv for input into geneious.")
            print("\n")
            if self.conf.print_result_to_screen == 1:
                input("press enter to continue.")
                self.intro_clear()
                print("Results will appear on screen as the file is processed.")
            else:
                self.intro_clear()
                print("Result are configured to only output to the results file. If would the the result to print to screen set \"print_result_to_screen\" in the dna.conf file, and restart the program.")

    def gen_id(self):
        _id = self.current_header.replace('@', '').split(':')
        return "{0}-{1}-{2}m".format(_id[0], _id[5], _id[6])

    def repeater(self, l, miner):
        repeats = []
        for _ in range(round(int(self.conf.read_average_length / len(l))), miner - 2, -1):
            if miner == 0:
                break
            else:
                repeats.append(re.compile(l * _))
        return repeats

    def permutation_generator(self):
        rainbow = []
        for _ in self.conf.permutation_elements:
            res = product("TGAC", repeat=int(_))
            for i in res:
                if len(i) == 4:
                    if i[0] + i[1] == i[2] + i[3]:
                        continue
                    else:
                        rainbow.append("".join(i))
                elif len(i) == 3:
                    if i[0] == i[1] == i[2]:
                        continue
                    else:
                        rainbow.append("".join(i))
                elif len(i) == 2:
                    if i[0] == i[1]:
                        continue
                    else:
                        rainbow.append("".join(i))
        self.rainbow = tuple(rainbow)

    def matrix(self):
        matches = []
        for _ in self.rainbow:
            l = _.strip();
            if len(l) == 4: matches.append(
                [re.compile(l), len(l), self.conf.min_adjacent_tetra, self.conf.read_average_length, self.repeater(l, self.conf.min_adjacent_tetra)])
            if len(l) == 3: matches.append(
                [re.compile(l), len(l), self.conf.min_adjacent_tri, self.conf.read_average_length, self.repeater(l, self.conf.min_adjacent_tri)])
            if len(l) == 2: matches.append(
                [re.compile(l), len(l), self.conf.min_adjacent_di, self.conf.read_average_length, self.repeater(l, self.conf.min_adjacent_di)])
        self.matches = matches

    def filter(self):
        for _ in self.matches:
            if self.current_dirty == 1:
                break
            elif _[1] == 4 and self.conf.filter_tetra_min_nonadjacent + self.conf.filter_tetra_min_adjacent > 0:
                if self.conf.filter_tetra_min_nonadjacent > 0 and len(re.findall(_[0], self.current_read_leadfilter)) > self.conf.filter_tetra_min_nonadjacent or len(
                        re.findall(_[0], self.current_read_rearfilter)) > self.conf.filter_tetra_min_nonadjacent:
                    self.current_dirty = 1
                    break
                elif self.conf.filter_tetra_min_adjacent > 0 and _[0].pattern * self.conf.filter_tetra_min_adjacent in self.current_read_leadfilter or _[
                    0].pattern * self.conf.filter_tetra_min_adjacent in self.current_read_rearfilter:
                    self.current_dirty = 1
                    break
            elif _[1] == 3 and self.conf.filter_tri_min_nonadjacent + self.conf.filter_tri_min_adjacent > 0:
                if self.conf.filter_tri_min_nonadjacent > 0 and len(re.findall(_[0], self.current_read_leadfilter)) > self.conf.filter_tri_min_nonadjacent or len(
                        re.findall(_[0], self.current_read_rearfilter)) > self.conf.filter_tri_min_nonadjacent:
                    self.current_dirty = 1
                    break
                elif self.conf.filter_tri_min_adjacent > 0 and _[0].pattern * self.conf.filter_tri_min_adjacent in self.current_read_leadfilter or _[
                    0].pattern * self.conf.filter_tri_min_adjacent in self.current_read_rearfilter:
                    self.current_dirty = 1
                    break
            elif _[1] == 2 and self.conf.filter_di_min_nonadjacent + self.conf.filter_di_min_adjacent > 0:
                if self.conf.filter_di_min_nonadjacent > 0 and len(re.findall(_[0], self.current_read_leadfilter)) > self.conf.filter_di_min_nonadjacent or len(
                        re.findall(_[0], self.current_read_rearfilter)) > self.conf.filter_di_min_nonadjacent:
                    self.current_dirty = 1
                    break
                elif self.conf.filter_di_min_adjacent > 0 and _[0].pattern * self.conf.filter_di_min_adjacent in self.current_read_leadfilter or _[
                    0].pattern * self.conf.filter_di_min_adjacent in self.current_read_rearfilter:
                    self.current_dirty = 1
                    break
            else:
                break

    def reporter(self):

        if self.current_id[-1:] == "p":
            report = "\nREAD ID: {0}\nFILE NAME: {1}\nLINE HDR: {2}\nN Count: {3}\nREAD Seq: {4}\n"
            return report.format(self.current_id, self.current_fn, self.current_header, len(re.findall("N", self.current_read)),
                             self.current_read)
        elif self.current_id[-1:] == "m":
            report = "\nREAD ID: {0}\nFILE NAME: {1}\nLINE HDR: {2}\nN Count: {3}\nResults: {4}\nREAD Seq: {5}\n"
            return report.format(self.current_id, self.current_fn, self.current_header, len(re.findall("N", self.current_read)),
                             self.current_results, self.current_read)

    def result_files(self):
        if self.current_id[-1:] == "m":
            with open(self.conf.result_csv, 'a') as csv_file:
                with open(self.conf.result_file, 'a') as result_file:
                    csv_file.write("{}, {}\n".format(self.current_header, self.current_read))
                    result_file.write(self.reporter())
        elif self.current_id[-1:] == "p":
            with open(self.conf.potential_results, 'a') as potential_file:
                potential_file.write(self.reporter())


def reader(argr):
    file_counter = 0
    for i in argr.conf.fn:
        if file_counter:
            argr.conf.output_files(num = 1)
        argr.current_counter = 0
        file_counter +=1
        with open(i, 'r') as file:
            for line in file:
                if '@HWI' in line:
                    argr.cleaner()
                    argr.current_fn = file.name
                    argr.current_header = line.replace(" ", "").rstrip('\n')
                    argr.current_read = file.readline().replace(" ", "").rstrip('\n')
                    argr.current_counter = argr.current_counter + 1
                    argr.current_dirty = 0
                    if argr.current_counter % 500000 == 0: print("{} lines have been processed.".format(argr.current_counter))
                    argr.current_id = argr.gen_id()
                    evaluation(argr)


def evaluation(argr):
    argr.filter()
    if argr.current_dirty == 0:
        for _ in argr.matches:
            flag = 0
            if argr.conf.return_di == 'no' and _[1] == 2: continue
            elif _[2] == 0: continue
            elif argr.conf.return_tri == 'no' and _[1] == 3: continue
            elif argr.conf.return_tetra == 'no' and _[1] == 4: continue
            argr.current_match_occur_in_read = len(re.findall(_[0], argr.current_read))
            argr.current_match_occur_outside_filter = len(re.findall(_[0], argr.current_read[argr.conf.lead_filter:-argr.conf.rear_filter]))
            if argr.current_match_occur_in_read > _[2] and not argr.current_dirty:
                for repeat in _[4]:
                    if argr.current_dirty == 1:
                        break
                    elif len(repeat.pattern) / _[1] > argr.current_match_occur_in_read:
                        continue
                    elif len(repeat.pattern) / _[1] < _[2]:
                        if argr.conf.return_potential == "yes" and _[1] == 4 and argr.current_match_occur_outside_filter > argr.conf.min_non_adjacent_tetra:
                            argr.current_dirty = 1
                            argr.current_id += "p"
                            if len(re.findall(_[0], argr.current_read[:argr.conf.lead_filter])) < argr.conf.max_non_adjacent_filter_tetra:
                                if len(re.findall(_[0], argr.current_read[-argr.conf.rear_filter:])) < argr.conf.max_non_adjacent_filter_tetra:
                                    argr.current_result = 1
                                    argr.current_pos_results.update({_[0].pattern: [argr.current_match_occur_outside_filter]})
                                    break
                                else:
                                    break
                            else:
                                break
                        elif re.search(repeat, argr.current_read):
                            argr.current_dirty = 1
                            for m in repeat.finditer(argr.current_read):
                                if m.span()[0] > argr.conf.lead_filter:
                                    if m.span()[1] < argr.conf.read_average_length - argr.conf.rear_filter:
                                        argr.current_result = 1
                                        argr.current_id += "m"
                                        argr.current_results.update({_[0].pattern: [argr.current_match_occur_in_read, len(repeat.pattern) / _[1]]})
                                    else:
                                        break
                                else:
                                    break
            else:
                continue
        if len(argr.current_results.keys()) > 0 and argr.current_result == 1:
            if argr.conf.print_result_to_screen == 1:
                print(argr.reporter())
            argr.result_files()
        elif len(argr.current_pos_results.keys()) > 0 and argr.current_result == 1:
            argr.result_files()

def main():
    print("\n")
    input("Ensure that the dna.py and dna.conf are in the same folder. \npress enter to continue\n\n\n\n")
    argr = Matcher(ConfRead())
    argr.intro(1)
    print(argr.conf.current_conf())
    reader(argr)


if __name__ == '__main__':
    main()
