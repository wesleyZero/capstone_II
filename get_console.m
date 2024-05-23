


function console = get_console() 

    console.divider = ...
    "_________________________________________________________________________";

    console.section = @console_section;
    console.subsection = @console_subsection;
end

function void = console_section(section_name)

    divider = "_________________________________________________________________________";

    fprintf("%s%s\n", section_name, divider);
end


function void = console_subsection(section_name, indent)

    divider = "_________________________________________________________________________";

    tabs = "";
    for i = 1:indent
        tabs = tabs + "\t";
    end
    
    fprintf(tabs + "%s%s\n",section_name, divider);
end