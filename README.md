# Capstone II
**UC Santa Barbara Chemical Engineering 2024**

**This readme is under construction 🚧**

### Table of Contents 

 * [Description of the Project](#Description)
 * [How do you simulate a chemical process?](#How_do_you_simulate_a_chemical_process?)
    * [Process Flow Diagram](#Process_Flow_Diagram)
    * [Reactor Model](#Reactor_model)
        * [Reaction Chemistry](#Reaction_Chemistry)
 * [Nomenclature](#Nomenclature)
    * [Supercritical Fluids](#Supercritical_fluids)
 * [Links](#Links)


### Description
This description and readme is a _massive_ simplification of our chemical engineering capstone for the intended audience of software engineers, with the entire emphasis being of how the code was written. For the full detailed report (with the intended audience of chemical engineers) please see the [PDF of our final Report](https://github.com/wesleyZero/capstone_II/blob/main/readme/pdf/capstone2_DMC.pdf).

We (my capstone group) performed a detailed techno-economic analysis and design for a proposed 100 kiloton per year (kta) plant producing polymer-grade Dimethyl Carbonate (DMC) for optical applications. Our economic analysis included (but is certainly not limited to) calculations of net present value (NPV), internal rate of return (IRR), total capital investment (TCI), energy consumption, carbon emissions, hazard ananylsis, and _much_ more. Those details will be ommitted in this readme since the focus of the following will be the code. 

# How do you simulate a chemical reactor?

You start off with a process flow diagram, this process flow diagram shows all processes that are included in the chemical processing plant that we are designing.

rough draft of logic flow


- overall idea
- process flow diagram (The big picture)
- zoom in of reactor (the focus)
- the primary function (root caller, call structure) 
- reaction chemistry
- chemical kinetics 
   - rate constant
   - concentration
      - flowrate
   - basis to plant scaling 



## Process Flow Diagram 

![Process Flow Diagram](https://github.com/wesleyZero/capstone_II/blob/main/readme/img/process_flow_diagram.jpeg)


**You can see above** all of the unit opererations (reactors, separators, distillation towers, etc.). **This readme will only be focusing on the reactor** near the top center-left of the diagram. **Why?** Well, in short chemical separations involve lots of "fudge factors" i.e. correlation factors and empirical models that account for the differences between simulation and reality. With a separation process this complex, it's not very useful to simulate the entire separation since the amount of error that will accumulate throughout the process is so large that it's neccessary to employ professional process simulation software like the one we used (Aspen HYSYS). Even the process simulation software has a large degree of error and requires a very educated chemical engineer to use correctly, and _even then_ there will be error.

## Reactor Model

We need a mathematical model for the reactor. The reactor we are using for this process is called a CSTR (constantly stirred tank reactor). It's kinda like if you imagine a witch brewing an evil potion or something like that 🧙‍♂️, however there are pipes that are flowing in and out of the brew constantly flowing reactants in and products out, there are lots of differential equations, and in this particular case some of the reactants are super critical fluids! Simple right?!

### Overview of how the simulation will work

**The end goal of this simulation is to answer just a few questions.** 
- What size of a reactor are we going to use?
- Are we going to operate the reactor isothermally? or isobarically? (i.e. at constant temp or pressure?)
     -If isothermal, what pressure are we going to operate at? (note: we were only allowed one isothermal temperature)
     -If isobaric, what temperature are we going to operate at? (note: we were only allowed one isobaric pressure to operate at)
- What is the economic value of the reactor that we designed? how do we choose conditions to optimize the economics of the process?

**How will we answer these questions?**
- Look at the reaction chemistry (what reactions can occur?, what species are involved?)
- Evaluate the chemical kinetics (given a set of conditions, how much of each reaction will occur?)
- 

**But Wes...What the heck is a super-critical fluid?**

Read in the Nomenclature section for [Supercritical Fluids](#supercritical_fluids) to learn what this is, if you don't know!

The following is going to be a bit hard to follow, whats important to know is **the end goal is to produce a set of graphs**. These graphs are going be be **functions of the residence time** (the mean time a chemical species spends in the reactor, symbol: tau), **conversion** (how much of the limiting species gets converted into products, symbol: chi), **temperature**, and **pressure**. The other important thing to note, is that we initially didn't know if an isothermal (constant temp) or isobaric (constant pressure) model was going to work. So we had to try out both, you will see this in the code.

## Reaction Chemistry

First, we need to know what chemical species are reacting, and how much. So we look at the reaction chemistry. 

![](https://github.com/wesleyZero/capstone_II/blob/main/readme/img/reaction_chemistry.png)

## Chemical Kinetics
Chemical kinetics, simply means "how fast is this reaction happening?" i.e. how much will the concentration of a species change per second?. To determine this, we must first look at the rate constant. 

### Rate constant

The **rate constant** is a constant that tells you how fast a particular reaction will happen, for a given concentration of reactants. The rate constant is a function of temperature, because at higher tempertures we get more chemical collisions and a higher concentration of species with enough energy to react upon a collision.

[get_isobaric_rate_constant](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L323)
```matlab
function k = get_isobaric_rate_constant(reaction, T)
    % input: 
    %   T [ C ]
    % output:
    %   k [mol / L s ]

    const = get_constants();
    thermo = const.thermo;
    T = const.units.temperature.c_to_k(T); % [ K ]

    switch reaction
        case '2f'
            k = 6.69 * 10^2 * exp(-37200 / (thermo.R * T));
        case '2r'
            k = 1.19 * 10^4 * exp(-53700 / (thermo.R * T));
        case '3'
            k = 1.89 * 10^6 * exp(-82400 / (thermo.R * T));
        otherwise
            k = NaN;
            disp("ERROR: get_isobaric_rate_constant(): invalid reaction option");
    end
end
```

[get_isothermal_rate_constant](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L346)
```matlab
function k = get_isothermal_rate_constant(reaction, T, P)
    % input:
    %   P [ bar ]
    % These functions are from the research paper
    rho = get_supercritical_c02_density(T, P, 'isothermal');
    switch reaction
        case '2f'
            if rho > 246.82 % [g / L]
                k = (2.486 * 10^(-2)) - (4.943 * (10^(-5)) * rho);
            else
                k = (1.362 * 10^(-2)) - (1.569 * (10^(-6)) * rho);
            end
        case '2r'
            k = 0.01486 * rho^(-0.873);
        case '3'
            k = 3.014 * (10^(-4)) * exp(-5.99 * (10^(-3)) * rho);
        otherwise
            k = NaN;
            disp("ERROR: get_isothermal_rate_constant(): invalid reaction option")
    end
end
```

### Reaction Rate

Once we have the rate constants, we can calculate the reaction rate. 

[get_reaction_rate](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L207)
```matlab
function r = get_reaction_rate(C, reaction, T, P, opt)
    % input:
    %   opt = 'isothermal' or 'isobaric'

    k = get_rate_constant(reaction, T, P, opt);
    switch reaction
        case '2f'
            r = k * (C.ethylene_carbonate)^0.8;
        case '2r'
            r = k * C.dimethyl_carbonate * C.ethylene_carbonate;
        case '3'
            r = k * C.ethylene_carbonate;
        otherwise
            r = NaN;
            disp("ERROR: get_reaction_rate(): invalid reaction option")
    end
end
```

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
