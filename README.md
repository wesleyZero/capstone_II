# Capstone II
**UC Santa Barbara Chemical Engineering 2024**

### Table of Contents 

 * [Description of the Project](#Description)
 * [How do you simulate a chemical process?](#How_do_you_simulate_a_chemical_process?)
    * [Process Flow Diagram](#Process_Flow_Diagram)
 * [Links](#Links)


### Description
This description and readme is a _massive_ simplification of our chemical engineering capstone for the intended audience of software engineers, with the entire emphasis being of how the code was written. For the full detailed report (with the intended audience of chemical engineers) please see the [PDF of our final Report](https://github.com/wesleyZero/capstone_II/blob/main/readme/pdf/capstone2_DMC.pdf).

We (my capstone group) performed a detailed techno-economic analysis and design for a proposed 100 kiloton per year (kta) plant producing polymer-grade Dimethyl Carbonate (DMC) for optical applications. Our economic analysis included (but is certainly not limited to) calculations of net present value (NPV), internal rate of return (IRR), total capital investment (TCI), energy consumption, carbon emissions, hazard ananylsis, and _much_ more. Those details will be ommitted in this readme since the focus of the following will be the code. 

## How_do_you_simulate_a_chemical_process?

You start off with a process flow diagram, this process flow diagram shows all processes that are included in the chemical processing plant that we are designing.

### Process_Flow_Diagram 

![Process Flow Diagram](https://github.com/wesleyZero/capstone_II/blob/main/readme/img/process_flow_diagram.jpeg)

**You can see above** all of the unit opererations (reactors, separators, distillation towers, etc.) This readme will only be focusing on the reactor near the top center-left of the diagram. Why? Well, in short chemical separations involve lots of "fudge factors" i.e. correlation factors and empirical models that account for the differences between simulation and reality. With a separation process this complex, it's not very useful to simulate the entire separation since the amount of error that will accumulate throughout the process is so large that it's neccessary to employ professional process simulation software like the one we used (Aspen HYSYS). Even the process simulation software has a large degree of error and requires a very educated chemical engineer to use correctly. 

### Reactor_Model

We need a mathematical model for the reactor. The reactor we are using for this process is called a CSTR (constantly stirred tank reactor). It's kinda like if you imagine a witch brewing an evil potion or something like that 🧙‍♂️, however there are pipes that are flowing in and out of the brew constantly flowing reactants in and products out, there are lots of differential equations and in this particular case some of the reactants are super critical fluids! Simple right?!



## Links
[PDF of our final Report](https://github.com/wesleyZero/capstone_II/blob/main/readme/pdf/capstone2_DMC.pdf)

