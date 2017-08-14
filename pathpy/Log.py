# -*- coding: utf-8 -*-
"""
    pathpy is an OpenSource python package for the analysis of sequential data
    on pathways and temporal networks using higher- and multi order graphical models

    Copyright (C) 2016-2017 Ingo Scholtes, ETH Zürich

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact the developer:

    E-mail: ischoltes@ethz.ch
    Web:    http://www.ingoscholtes.net
"""

import enum
import time
import io
import sys

class Severity(enum.IntEnum):
    """ An enumeration that can be used to indicate
        the severity of log messages, and which can be
        used tpo filter messages based on severities.
    """

    ## Error messages
    ERROR = 4

    ## Warning messages
    WARNING = 3

    ## Informational messages (default minimum level)
    INFO = 2

    ## Messages regarding timing and performance
    TIMING = 1

    ## Debug messages (really verbose)
    DEBUG = 0


class Log:
    """ A simple logging class, that allows to select what messages should
        be recorded in the output, and where these message should be directed.
    """

    ## the output stream to which log entries will be written
    output_stream = sys.stdout

    ## The minimum severity level of messages to be logged
    min_severity  = Severity.INFO


    @staticmethod
    def setMinSeverity(severity):
        """ Sets the minimum sveerity level a message
        needs to have in order to be recorded in the output stream.
        By default, any message which has a severity of at least
        Severity.INFO will be written to the output stream. All messages
        with lower priority will be surpressed.
        """
        Log.min_severity = severity


    @staticmethod
    def setOutputStream(stream):
        """ Sets the output stream to which all messages will be
            written. By default, this is sys.stdout, but it can be
            changed in order to redirect the log to a logfile.
        """
        output_stream = stream


    @staticmethod
    def add(msg, severity=Severity.INFO):
        """ Adds a message with the given severity to the log. This message will be written
            to the log output stream, which by default is sys.stdout. A newline character
            will be added to the message by default.
        """
        if severity >= Log.min_severity:
            ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            Log.output_stream.write(ts + ' [' + str(severity) + ']\t' + msg + '\n')
            Log.output_stream.flush()
