# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Massimiliano Volpe & Lauri Mesilaakso"
__copyright__ = "Copyright 2022, Massimiliano Volpe & Lauri Mesilaakso"
__email__ = "massimiliano.volpe@liu.se & lauri.mesilaakso@regionostergotland.se"
__license__ = "GPL-3"


def test_dummy():
    from dummy import dummy
    assert dummy() == 1
