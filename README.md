# Capstone II
**UC Santa Barbara Chemical Engineering 2024**

### Table of Contents 

 * [Description of the Project](#Description)
 * [How do you simulate a chemical process?](#How_do_you_simulate_a_chemical_process?)
    * [Process Flow Diagram](#Process_Flow_Diagram)
    * [Reactor Model](#Reactor_model)
 * [Nomenclature](#Nomenclature)
    * [Supercritical Fluids](#Supercritical_fluids)
 * [Links](#Links)


### Description
This description and readme is a _massive_ simplification of our chemical engineering capstone for the intended audience of software engineers, with the entire emphasis being of how the code was written. For the full detailed report (with the intended audience of chemical engineers) please see the [PDF of our final Report](https://github.com/wesleyZero/capstone_II/blob/main/readme/pdf/capstone2_DMC.pdf).

We (my capstone group) performed a detailed techno-economic analysis and design for a proposed 100 kiloton per year (kta) plant producing polymer-grade Dimethyl Carbonate (DMC) for optical applications. Our economic analysis included (but is certainly not limited to) calculations of net present value (NPV), internal rate of return (IRR), total capital investment (TCI), energy consumption, carbon emissions, hazard ananylsis, and _much_ more. Those details will be ommitted in this readme since the focus of the following will be the code. 

# How_do_you_simulate_a_chemical_process?

You start off with a process flow diagram, this process flow diagram shows all processes that are included in the chemical processing plant that we are designing.

# Process_Flow_Diagram 

![Process Flow Diagram](https://github.com/wesleyZero/capstone_II/blob/main/readme/img/process_flow_diagram.jpeg)


**You can see above** all of the unit opererations (reactors, separators, distillation towers, etc.) This readme will only be focusing on the reactor near the top center-left of the diagram. Why? Well, in short chemical separations involve lots of "fudge factors" i.e. correlation factors and empirical models that account for the differences between simulation and reality. With a separation process this complex, it's not very useful to simulate the entire separation since the amount of error that will accumulate throughout the process is so large that it's neccessary to employ professional process simulation software like the one we used (Aspen HYSYS). Even the process simulation software has a large degree of error and requires a very educated chemical engineer to use correctly, and _even then_ there will be error.

# Reactor_Model

We need a mathematical model for the reactor. The reactor we are using for this process is called a CSTR (constantly stirred tank reactor). It's kinda like if you imagine a witch brewing an evil potion or something like that 🧙‍♂️, however there are pipes that are flowing in and out of the brew constantly flowing reactants in and products out, there are lots of differential equations and in this particular case some of the reactants are super critical fluids! Simple right?!

**But Wes...What the heck is a super-critical fluid?**

Read in the Nomenclature section for [Supercritical Fluids](#supercritical_fluids) to learn what this is, if you don't know!

The following is going to be a bit hard to follow, whats important to know is **the end goal is to produce a set of graphs**. These graphs are going be be **functions of the residence time** (the mean time a chemical species spends in the reactor, symbol: tau), **conversion** (how much of the limiting species gets converted into products, symbol: chi), **temperature**, and **pressure**. 

# Nomenclature

## Supercritical_fluids
<img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/supercritical_phase_diagram.png" width="500">

**A supercritical fluid** is a fluid that behaves like a gas and a liquid simultaneously. Depending on the thermodynamic conditions (the temperature and pressure) it will behave more like one or the other. The details of this are not what I want the focus to be, thats a topic of physics. Whats important is that we have a reactive component that changes density, **effectively changing the concentration of all of the species in the reactor depending on what thermodynamic conditions we specifiy**. For CO2 we get this type of bizarro phase behavior when we get to pressures of about 76 times atmospheric. 

The function to model this is below. The functions you see are empirical regressions of experimental data since there generally isn't good physical models of these exotic phases of matter. 

[get_supercritical_c02_density function](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L368)
```matlab
function rho = get_supercritical_c02_density(T, P, opt)
    % Input: condition = T or P. Depending on option
    %   P [=] bar
    %   T [=] celcius   
    % Assumptions:
    %   Isobaric model is at 150 bar
    %   Isothermal model is at 140 C
    % Ranges of input
    %   P = [50 bar, 150 bar]
    %   T = [80 C, 140 C]
    % Output: 
    %   rho [=] kg / m^3
    withinTempRange = @(T) T >= 80 && T <= 140;
    withinPressureRange = @(P) P >= 50 && P <= 150;
    
    rho.units = 'kg / m^3';
    switch opt
        case 'isothermal'
            if withinPressureRange(P)
                rho = 1.6746 * P - 12.592;
                    % NIST / Excel Regression
            else
                rho = NaN;
                disp("get_supercritical_c02_density : ERROR : P out of range")
            end
        case 'isobaric'
            if withinTempRange(T)
                if P < 125 % [ Bar ]
                    rho = 356.08 * exp(-0.006 * T);
                        % NIST Data at 100 bar
                else
                    rho = 838.87 * exp(-0.009 * T);
                        % NIST Data at 150 bar
                end
            else
                rho = NaN;
                disp("get_supercritical_c02_density : ERROR : T out of range")
            end
        otherwise
            disp("SUPERCRITICAL C02 DENSITY FUNCTION ERROR: invalid opt")
            rho = NaN;
    end
end
```

## Links
[PDF of our final Report](https://github.com/wesleyZero/capstone_II/blob/main/readme/pdf/capstone2_DMC.pdf)
