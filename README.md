# MRCatheterTrackingSim
Tracking Simulator Client



# reformat_log.py

This is a script to reformat a log file. Due to a bug in the sequence,
the 11th and 12 th numbers are concatenated. The script will split them
assuming that the number of digits in each column is six. Note that
the number of digits may be smaller if the last digit is zero, and
this assumption may not always be true.

