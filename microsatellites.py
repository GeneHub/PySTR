#!/usr/bin/python3.4

######################################################
#Author : Robert Battocletti
#Program Name: microsatellite finder
#Date Published: 26FEB2014
#Designed to evaluate fastq files for microsatellites
#This update includes multithreading
######################################################
import re
from itertools import product
import configparser
from glob import glob
from sys import exit
import threading
from time import localtime as time


#Globals
line_num = 0
mfile = ""
exitFlag = 0
threadLock0 = threading.Lock()
threadLock1 = threading.Lock()


class myThread (threading.Thread):
    def __init__(self, threadID, name, counter,argr):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.counter = counter
        self.argr = argr
        self.argr.current_read = "primer"
    def run(self):
        print("Starting " + self.name)
        while not exitFlag:
            threadLock0.acquire()
            reader(self.argr)
            threadLock0.release()
            if self.argr.current_read and not exitFlag:
                self.run_filter_test()
            else:
                break
        print("Exiting " + self.name)
    def run_filter_test(self):
        self.argr.filter()
        if self.argr.current_dirty == 0:
            self.test_evaluation()
    def test_evaluation(self):
        evaluation(self.argr)
        if self.argr.current_result == 1 and self.argr.current_dirty == 0:
            self.result_call()
    def result_call(self):
        threadLock1.acquire()
        self.argr.result_files()
        threadLock1.release()


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
        #self.project_name = config.get('user_defined', 'project_name')
        self.project_name = ""
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
        #self.return_potential = str(config.get('user_defined', 'return_potential')).lower() #Disabled requires filter work.
        self.return_potential = "no"
        self.process_directory = int(config.get('computing', 'process_directory'))
        self.print_result_to_screen = int(config.get('computing', 'print_result_to_screen'))
        self.max_non_adjacent_filter_tetra = int(config.get('user_defined', 'max_non_adjacent_filter_tetra'))
        self.read_average_length = int(151)
        self.write_database = "no"
        self.result_file = self.project_name + "_microsatellites.txt"
        self.database = self.project_name + "_.db"
        self.result_csv = self.project_name + "_results.csv"
        self.potential_results = self.project_name + "_potential.txt"
        self.fn = []

    def __repr__(self):
        return repr((self.project_name, self.fn))

    def current_conf(self):
        cconf = '\n----------------\nCurrent configuration summary:\n----------------\nProject Name: {0}\nProcess folder:{7}\n--------------\nReturn tetra: {1}\nMin_adjacent tetra: {2}\nReturn tri: {3}\nMin_adjacent tri: {4}\nReturn di: {5}\nMin adjacent di: {6}\n--------------\nLead filter area: {8}\nRear filter area: {9}\nFilter tetra min non-adjacent: {10}\nFilter tetra min adjacent: {11}\nFilter tri min non-adjacent: {12}\nfilter tri min adjacent: {13}\nfilter di min non-adjacent:{14}\nfilter di min adjacent: {15}\n----------------\n'
        return cconf.format(self.project_name, self.return_tetra, self.min_adjacent_tetra, self.return_tri,
                            self.min_adjacent_tri, self.return_di, self.min_adjacent_di, self.process_directory,
                            self.lead_filter,
                            self.rear_filter, self.filter_tetra_min_nonadjacent, self.filter_tetra_min_adjacent,
                            self.filter_tri_min_nonadjacent, self.filter_tri_min_adjacent,
                            self.filter_di_min_nonadjacent, self.filter_di_min_adjacent)

    def output_files(self, num=0):
        try:
            with open(self.project_name + self.result_file, "a") as result_file:
                result_file.write(self.current_conf())
        except Exception as exception:
            print(exception)
        if num == 0:
            if self.write_database == "yes":
                try:
                    with open(self.project_name +self.database, "a"):
                        pass
                except Exception as exception:
                    print(exception)
            try:
                with open(self.project_name +self.result_csv, "a") as csv_file:
                    csv_file.write("Name,Sequence\n")
            except Exception as exception:
                print(exception)
            if self.return_potential == "yes":
                try:
                    with open(self.project_name +self.potential_results, "a"):
                        pass
                except Exception as exception:
                    print(exception)

    def read_directory(self):
        fn = []
        if self.process_directory == 0:
            fn.append(input("Enter the file you would like to process :> "))
            #fn.append("lane2_Undetermined_L002_R1_033.fastq")
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
                    fn = fn_tmp;
                    switch = 1
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
        self.current_results = {}
        self.current_pos_results = {}
        self.current_counter = int(0)
        self.current_result = 0
        self.current_dirty = 0
        self.current_match_occur_in_read = 0
        self.current_match_occur_outside_filter = 0

    def filter_define(self):
        self.current_read_leadfilter = self.current_read[:self.conf.lead_filter]
        self.current_read_rearfilter = self.current_read[-self.conf.rear_filter:]

    def cleaner(self):
        self.current_fn = ''
        self.current_header = ''
        self.current_read = ''
        self.current_id = ''
        self.current_results = {}
        self.current_pos_results = {}
        self.current_result = 0
        self.current_dirty = 0
        self.current_pos_result = 0

    def intro(self, num):
        if num == 0:
            print("\n\n")
            input("Ensure that the dna.py and dna.conf are in the same folder. \npress enter to continue\n\n")
        elif num == 1:
            print("This process could take up to an one hour per 16GB file.")
            print("When complete, there will three different output files:")
            print("one database, one text file with results and a csv for input into geneious.")
            print("\n")
            if self.conf.print_result_to_screen == 1:
                input("press enter to continue.")
                print("Results will appear on screen as the file is processed.")
            else:
                print(
                    "Result are configured to only output to the results file. If would the the result to print to screen set \"print_result_to_screen\" in the dna.conf file, and restart the program.")

    def gen_id(self):
        _id = self.current_header.replace('@', '').split(':')
        return "{0}-{1}-{2}".format(_id[0], _id[5], _id[6])

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
                rainbow.append("".join(i))
        self.rainbow = tuple(rainbow)

    def matrix(self):
        matches = []
        for _ in self.rainbow:
            filter_only = False
            l = _.strip()
            if len(l) == 4:
                if _[0] + _[1] == _[2] + _[3]:
                    filter_only = True
                matches.append([re.compile(l), len(l), self.conf.min_adjacent_tetra, self.conf.read_average_length,
                                self.repeater(l, self.conf.min_adjacent_tetra), filter_only])
            if len(l) == 3:
                if _[0] == _[1] == _[2]:
                    filter_only = True
                matches.append([re.compile(l), len(l), self.conf.min_adjacent_tri, self.conf.read_average_length,
                                self.repeater(l, self.conf.min_adjacent_tri), filter_only])
            if len(l) == 2:
                if _[0] == _[1]:
                    filter_only = True
                matches.append([re.compile(l), len(l), self.conf.min_adjacent_di, self.conf.read_average_length,
                                self.repeater(l, self.conf.min_adjacent_di), filter_only])
        self.matches = matches

    def filter(self):
        self.filter_define()
        for _ in self.matches:
            if self.current_dirty == 1:
                break
            elif _[1] == 4 and self.conf.filter_tetra_min_nonadjacent + self.conf.filter_tetra_min_adjacent > 0:
                if self.conf.filter_tetra_min_nonadjacent > 0 and len(
                        re.findall(_[0], self.current_read_leadfilter)) > self.conf.filter_tetra_min_nonadjacent or len(
                        re.findall(_[0], self.current_read_rearfilter)) > self.conf.filter_tetra_min_nonadjacent:
                    self.current_dirty = 1
                    break
                elif self.conf.filter_tetra_min_adjacent > 0 and _[
                    0].pattern * self.conf.filter_tetra_min_adjacent in self.current_read_leadfilter or _[
                    0].pattern * self.conf.filter_tetra_min_adjacent in self.current_read_rearfilter:
                    self.current_dirty = 1
                    break

            elif _[1] == 3 and self.conf.filter_tri_min_nonadjacent + self.conf.filter_tri_min_adjacent > 0:
                if self.conf.filter_tri_min_nonadjacent > 0 and len(
                        re.findall(_[0], self.current_read_leadfilter)) > self.conf.filter_tri_min_nonadjacent or len(
                        re.findall(_[0], self.current_read_rearfilter)) > self.conf.filter_tri_min_nonadjacent:
                    self.current_dirty = 1
                    break
                elif self.conf.filter_tri_min_adjacent > 0 and _[
                    0].pattern * self.conf.filter_tri_min_adjacent in self.current_read_leadfilter or _[
                    0].pattern * self.conf.filter_tri_min_adjacent in self.current_read_rearfilter:
                    self.current_dirty = 1
                    break

            elif _[1] == 2 and self.conf.filter_di_min_nonadjacent + self.conf.filter_di_min_adjacent > 0:
                if self.conf.filter_di_min_nonadjacent > 0 and len(
                        re.findall(_[0], self.current_read_leadfilter)) > self.conf.filter_di_min_nonadjacent or len(
                        re.findall(_[0], self.current_read_rearfilter)) > self.conf.filter_di_min_nonadjacent:
                    self.current_dirty = 1
                    break
                elif self.conf.filter_di_min_adjacent > 0 and _[
                    0].pattern * self.conf.filter_di_min_adjacent in self.current_read_leadfilter or _[
                    0].pattern * self.conf.filter_di_min_adjacent in self.current_read_rearfilter:
                    self.current_dirty = 1
                    break

    def reporter(self):
        if self.current_result == 0:
            report = "\nREAD ID: {0}\nFILE NAME: {1}\nLINE HDR: {2}\nN Count: {3}\nREAD Seq: {4}\n"
            return report.format(self.current_id, self.current_fn, self.current_header,
                                 len(re.findall("N", self.current_read)),
                                 self.current_read)
        elif self.current_result > 0:
            report = "\nREAD ID: {0}\nFILE NAME: {1}\nLINE HDR: {2}\nN Count: {3}\nResults: {4}\nREAD Seq: {5}\n"
            return report.format(self.current_id, self.current_fn, self.current_header,
                                 len(re.findall("N", self.current_read)),
                                 self.current_results, self.current_read)

    def result_files(self):
        if self.current_result > 0:
            with open(self.conf.project_name + self.conf.result_csv, 'a') as csv_file:
                with open(self.conf.project_name + self.conf.result_file, 'a') as result_file:
                    csv_file.write("{}, {}\n".format(self.current_header, self.current_read))
                    result_file.write(self.reporter())
        elif self.current_result == 0:
            with open(self.conf.project_name + self.conf.potential_results, 'a') as potential_file:
                potential_file.write(self.reporter())

def reader(argr):
    global mfile, line_num, exitFlag
    read = True
    while read:
        global mfile
        tmp_1 = next(mfile, '')
        if tmp_1:
            if tmp_1[:4] == '@HWI':
                argr.cleaner()
                argr.current_fn = mfile.name
                argr.current_header = tmp_1.replace(" ", "").rstrip('\n')
                argr.current_read = next(mfile).replace(" ", "").rstrip('\n')
                line_num += 2
                argr.current_counter += 1
                if line_num % 200000 == 0:
                    progress_time = [ "at ", time()]
                    print("{}% has been processed.{}".format(round(( line_num / 4000000) *100)), print_time(progress_time))
                if len(re.findall("N", argr.current_read)) > 18:
                    continue
                argr.current_id = argr.gen_id()
                read = False
                break
        else:
            read = False
            exitFlag = 1
            break

def evaluation(argr):
    if argr.current_dirty == 0:
        for _ in argr.matches:
            if argr.conf.return_di == 'no' and _[1] == 2:
                continue
            elif _[5] == True:
                continue
            elif argr.conf.return_tri == 'no' and _[1] == 3:
                continue
            elif argr.conf.return_tetra == 'no' and _[1] == 4:
                continue
            elif _[2] == 0:
                continue
            argr.current_match_occur_in_read = len(re.findall(_[0], argr.current_read))
            argr.current_match_occur_outside_filter = len(
                re.findall(_[0], argr.current_read[argr.conf.lead_filter:-argr.conf.rear_filter]))
            if argr.current_match_occur_in_read > _[2] and not argr.current_dirty:
                for repeat in _[4]:
                    if argr.current_dirty == 1:
                        break
                    elif len(repeat.pattern) / _[1] > argr.current_match_occur_in_read:
                        continue
                    elif len(repeat.pattern) / _[1] < _[2]:
                        if argr.conf.return_potential == "yes" and _[
                            1] == 4 and argr.current_match_occur_outside_filter > argr.conf.min_non_adjacent_tetra:
                            if len(re.findall(_[0], argr.current_read[
                                                    :argr.conf.lead_filter + 5])) < argr.conf.max_non_adjacent_filter_tetra - 1:
                                if len(re.findall(_[0], argr.current_read[-(
                                            argr.conf.rear_filter + 5):])) < argr.conf.max_non_adjacent_filter_tetra - 1:
                                    argr.current_pos_result = 1
                                    argr.current_pos_results.update(
                                        {_[0].pattern: [argr.current_match_occur_outside_filter]})
                                    break
                                else:
                                    argr.current_dirty = 1
                                    break
                            else:
                                argr.current_dirty = 1
                                break
                        elif re.search(repeat, argr.current_read):
                            for m in repeat.finditer(argr.current_read):
                                if m.span()[0] > argr.conf.lead_filter:
                                    if m.span()[1] < argr.conf.read_average_length - argr.conf.rear_filter:
                                        argr.current_result = 1
                                        argr.current_results.update({
                                            _[0].pattern: [argr.current_match_occur_in_read, len(repeat.pattern) / _[1]]})
                                    else:
                                        break
                                else:
                                    break


            else:
                continue
        if argr.current_result == 1 and argr.current_dirty == 0:
            argr.current_id += "m"
            if argr.conf.print_result_to_screen == 1:
                print(argr.reporter())
        elif argr.conf.return_potential == "yes" and argr.current_pos_result and argr.current_dirty == 0:
            argr.current_id += "p"
            argr.result_files()
        else:
            argr.current_dirty = 1

def update_argrs(argrs, threads):
    argrs[0].conf.read_directory()
    for i in range(0, threads):
        argrs[i].conf.fn = argrs[0].conf.fn

def print_time(time_stamp):
    timer = "\n{} = {}:{}:{}\n"
    return timer.format(time_stamp[0], time_stamp[1].tm_hour, time_stamp[1].tm_min, time_stamp[1].tm_sec)

def main():
    argrs = []
    threadlist = []
    print("\n")
    input("Ensure that the microsatellites.py and microsatellites.conf are in the same folder. \npress enter to continue\n\n\n\n")
    print("Enter the number of threads you would like to use.")
    print("Note: it will take a couple seconds to initialize the threads.")
    print("Default is 4.")
    user_input = input(":>")
    if len(user_input) == 0:
        threads = 4
    else:
        threads = int(user_input)
    for i in range(threads):
        argr = Matcher(ConfRead())
        argrs.append(argr)
        thid = myThread(i, "thread_" + str(i), i, argr)
        threadlist.append(thid)
    update_argrs(argrs, threads)
    start_time = ["Start Time:", time()]
    print(argrs[0].conf.current_conf())
    print(print_time(start_time))
    for file in argrs[0].conf.fn:
        global mfile
        for arid in argrs:
            arid.conf.project_name = file
        argrs[0].conf.output_files()
        mfile = open(file, 'r')
        for thread in threadlist:
            thread.start()
        for thread in threadlist:
            thread.join()
    endtime = ["Stop Time:", time()]
    print(print_time(start_time))
    print(print_time(endtime))

if __name__ == '__main__':
    main()
