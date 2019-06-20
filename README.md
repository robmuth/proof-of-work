# P💥W 
A zero-knowledge Proof-of-Work (PoW) to make sure that working hours of an employee is conform to the European law, withtout knowing the details such as starting and finishing times or breaks.

## Background
To ensure "health and safety" employees in Europe, the European Court of Justice decided in a ruling ([Case number C-55/18](http://curia.europa.eu/juris/liste.jsf?num=C-55/18&language=EN)) that its mandatory for employers to track the working hours of its employees.

Although this ruling has been spoken in favor of employees, ironically it can be also be used by employers against them by profiling working hours of individuals and imposing sanctions which might be contrary to the goal of the ruling.

## Workflow
Three entities are required to realize our scheme:

  1. An employee, providing the zero-knowledge proof
  2. An emploter, verifying the proof
  3. A terminal issuing timestamps, trusted by both other parties

