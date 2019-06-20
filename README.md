# P💥W
A zero-knowledge Proof-of-Work (PoW) to make sure that working hours of an employee conform to the European law, without knowing the details such as starting and finishing times or breaks.

## Background
To ensure "health and safety" employees in Europe, the European Court of Justice decided in a ruling ([Case number C-55/18](http://curia.europa.eu/juris/liste.jsf?num=C-55/18&language=EN)) that its mandatory for employers to track the working hours of its employees.

Although this ruling has been spoken in favor of employees, ironically it can be also be used by employers against them by profiling working hours of individuals and imposing sanctions which might be contrary to the goal of the ruling.

## Requirements
Bash scripts are optimized for Unix and the behavior is unknown for Linux! For python scripts (e.g., [`terminal.py`](./blob/master/bin/terminal.py)), Python3 is required and [`pycrypto`](https://github.com/Zokrates/pycrypto) (of Zokrates project) must be on `PYTHONPATH`.

## Workflow

![Workflow](./workflow.png)

## Quickstart
