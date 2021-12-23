# Another related proyect
#  - https://scimath.readthedocs.io/en/latest/index.html

metricu = "kg,m,s,C,mol,K".split(",")
metric_order = {i:j for i,j in  zip( metricu, range(len(metricu)) ) }

# equivalences to SI/metric, e.g.  { 'nm' : 10e-9 (meteres) , 'mi' : 1.6e3 (meters) , ... }
convertions = {
    "s"  : {            # Second
    "s"  : 1,
    "year":3.15576e7,   # Year
    "day": 86400,       # Day    
    "hr" : 3600,        # Hour
    "min": 60,          # Minute
    "ms" : 1e-6,        # Milisec
    "us" : 1e-6,        # Microsec
    "ns" : 1e-9,        # Nanosec
    },
    "kg" : {            # Kilogram
    "kg" : 1,
    "ug" : 1e-9,        # Miligran
    "mg" : 1e-6,        # Miligran
    "g"  : 1e-3,        # gram
    "ton": 9.071847e2,  # ton
    "lb" : 0.4535924,   # pound
    },
    "m"  : {            # Meter
    "m"  : 1,
    "A"  : 1e-10,       # Angstrom
    "nm" : 1e-9,        # Nanometer
    "um" : 1e-6,        # Micrometer
    "mm" : 1e-3,        # Milimeter
    "cm" : 1e-2,        # Centemiter
    "km" : 1e3,         # Kilometer
    "in" : 0.0254,      # Inches
    "ft" : 0.3048,      # Feet
    "yd" : 0.9144,      # Yard
    "mi" : 1.609344e3,  # Mile
    "au" : 1.495979e11  # Astronomical unit
    },
    "mol" : {
    "mol" : 1,
    "particles" : 1/(6.02214076e23) # Particulas
    },
    "C" : {             #Coulomb
        "C" : 1
    },
    "K" : {             #Kelvin
        "K" : 1
    }
}

combined_units = {
    "V"  : "kg.m^2.s^-2.C^-1",    # Volts
    "Amp": "s^-1.C",              # Ampere
    "Ohm": "kg.m^2.s^-1.C^-2",    # Ohm
    "N"  : "kg.m.s^-2",           # Newton
    "J"  : "kg.m^2.s^-2",         # Joule
    "Pa" : "kg.m^-1.s^-2",        # Pascal
    "L"  : "m^3",                 # Liters
    "W"  : "kg.m^2.s^-3"          # Watts
}

mconv = {}  # this smt like {"g":"kg", "lb":"kg" ...} it maps [unit] -> [metric_unit]
for u_metric in convertions.keys(): #{"s","kg","m","mol", "C"}
    for i in convertions[u_metric].keys():
        mconv.update( {i :  u_metric} ) 

iconv = {}  # like mconv but it maps [unit] -> [imperial_unit] 
for i,u_imperial in enumerate(["s","lb","ft","mol","C"]): # {"s","lb","ft","mol","C" }
    u_metric = list(convertions.keys())
    for i in convertions[u_metric[i]].keys():
        iconv.update( {i :  u_imperial} ) 

# Main Class --------------------------------------------------------------------
class Unit(object):
    """
    Abstraction of a Unit, includes convertions and basic operations. \n
    ------ EXAMPLE : Free Fall (No air resistance) -------

    g = Constants.g \n
    velocity = Unit(1, "km.s^-1") \n
    time = Unit(2, "min") \n
    distance = - g*time**2*(1/2) + velocity * time * 2  

    print( distance )                   # Unit( 4.93680e+04 m )  \n
    print( distance.to("mi").text() )   # 30.67585 km 

    --------------------------------------------------------
    Notes:
     - When having a unit to the zero e.g. km^0 in an operation, it'll give a float, ussually
     - Cannot multiply by const on the left i.e. 2*Unit() gives Error, but Unit()*2 is valid
     - Try to use parenthesis, i.e. Unit()* 1/2 gives Error, but Unit* (1/2) is valid
     - Every operation translates back to metric, u have to convert back after it 

    Limitations:
     - Can't compute sqrt, sin, cos, exp, or Unit()**Unit(), cuz is not implemented
     - TODO : implement temperatures units (°C,°F,K,...)
     - TODO : add error computations
    """
    def __init__(self, value:float, unit_str:str):
        # Parsing the units ---------------------------
        temp = unit_str
        for cu in combined_units.keys(): # Check if it has "cu"'s : {J, N, pa, ...}
            if cu in unit_str:
                cu_index = unit_str.rfind(cu)         # Find the position of the "cu" 
                if unit_str.rfind(cu+"^") != -1:      # If has a power
                    f_index = unit_str[cu_index+len(cu):].rfind(".")
                    if f_index != -1:
                        exponent = int( unit_str[cu_index+len(cu)+1:cu_index+len(cu)+f_index] ) # get the power
                    else:
                        exponent = int( unit_str[cu_index+len(cu)+1:] )          # get the power
                    cu_list = self.get_u_matrix(combined_units[cu])
                    cu_list = [ [i[0],exponent*i[1]] for i in cu_list]
                    temp = temp.replace(cu+f"^{exponent}", self._u_text2(cu_list))
                else:
                    temp = temp.replace(cu,combined_units[cu])
        self._units = self.get_u_matrix(temp)
        # self._units has the form [["kg",1],["hr",4],...] 
        # may be repeated, cuz of <<Combined Units>> so
        self._value = value 
        self._s_simplify()      # Eliminate repeated values 
        self.to_metric()        # Convert to metric
        self._units.sort(key = lambda s:metric_order[s[0]]) # Just for kiks

    ### Utility functions
    def get_u_matrix(self, u_str:str) -> list:
        """Input:  "kg^2.s^-3.cm^5" | Ouput: [['kg',2],['s',-3],['cm',5]]"""
        res = u_str.split(".")
        for i,u in enumerate(res): 
            if "^" in u:
                index = u.rfind("^")
                exponent = int( u[index+1:] )
                res[i] = [u[:index], exponent]
            else:
                res[i] = [u, 1]
        return res
    def get_units(self) -> list:
        return self._units
    def get_value(self) -> float:
        return self._value
    def _u_text(self) -> str:
        u_text = ""
        if self._units == None:
            return u_text
        for uni in self._units:
            if uni[1] != 1:
                u_text += f"{uni[0]}^{uni[1]}."
            else:
                u_text += f"{uni[0]}."
        return u_text[:-1]
    def _u_text2(self, u_list:list) -> str:
            u_text = ""
            if u_list == None:
                return u_text
            for uni in u_list:
                if uni[1] != 1:
                    u_text += f"{uni[0]}^{uni[1]}."
                else:
                    u_text += f"{uni[0]}."
            return u_text[:-1]
    def _simplify(self, units_l:list):
        temp = set( (i[0] for i in units_l) )
        u_list = [[i,0] for i in temp]
        for u in units_l:
            for i,un in enumerate(temp):
                if (u[0] == un):
                    u_list[i][1] += u[1]
                    break
        u_list = [i for i in u_list if i[1] != 0]
        if u_list != None:
            u_list.sort(key = lambda s:metric_order[s[0]])
        return u_list
    def _s_simplify(self):
        temp = set( (i[0] for i in self.get_units()) )
        u_list = [[i,0] for i in temp]
        for u in self.get_units() :
            for i,un in enumerate(temp):
                if (u[0] == un):
                    u_list[i][1] += u[1]
                    break
        u_list = [i for i in u_list if i[1] != 0]
        self._units = u_list
    def to_metric(self): # Shorthand 
        self.to("metric")
    def to_imperial(self): # Shorthand 
        return self.to("imperial")
    def n(self): # Shorthand
        return self._value
    def __repr__(self) -> str:
        return "Unit( {:.5e} {:s} )".format(self._value, self._u_text())
    def text(self, style="normal") -> str:
        if style == "sci":
            return "{:.8e} {:s}".format(self._value, self._u_text())
        elif style == "normal":
            return "{:.5f} {}".format(self._value, self._u_text())
        else:
            return "{:.9e} {:s}".format(self._value, self._u_text())
    def has_mathching_units(self, other) -> bool:
        other.to_metric()
        self.to_metric()
        return (self._u_text() == other._u_text())

    ### Mathematical Operations
    def __add__(self, other):
        if self.has_mathching_units(other):
            return Unit(self._value + other._value, self._u_text() )
        else:
            raise Exception(f"Can't add {self._u_text()} to {other._u_text()}")
    def __sub__(self, other):
        other.to_metric()
        self.to_metric()
        if (self._u_text() == other._u_text()):
            return Unit(self._value - other._value, self._u_text() )
        else:
            raise Exception(f"Can't subtract {self._u_text()} to {other._u_text()}")
    def __mul__(self, other):
        self.to_metric()
        if (type(other) == float) or (type(other) == int):
            return Unit(self._value * other, self._u_text() )
        else:
            other.to_metric()
            u_sum = self._units + other._units
            # print(u_sum)
            all_units = self._simplify(self._units + other._units)
            u_text = self._u_text2(all_units)
            if(u_text == ""):
                return self._value * other._value
            else:
                return Unit(self._value * other._value, u_text )
    def __truediv__(self, other):
        #inverting the units
        u_list = other._units
        for u in u_list:
            u[1] *= -1
        #using multiply to handle the stuff
        temp = Unit(1/other._value, self._u_text2(u_list))
        return self * temp
    def __lt__(self, other) -> bool:
        if self.has_mathching_units(other):
            return (self.get_value() < other.get_value())
        else:
            raise Exception(f"Can't compare {self._u_text()} to {other._u_text()}")
    def __gt__(self, other) -> bool:
        if self.has_mathching_units(other):
            return (self.get_value() > other.get_value())
        else:
            raise Exception(f"Can't compare {self._u_text()} to {other._u_text()}")
    def __eq__(self, other) -> bool: # Aproximately equal, that is
        if(self.has_mathching_units(other)):
            a = self.get_value()
            b = other.get_value()
            return (abs(a-b)/max(a,b) < 1e-10)
        else:
            raise Exception(f"Can't compare {self._u_text()} to {other._u_text()}")
    def compare(self, other) -> float:
        """
        If possible, computes the value of this/other. 
        e.g. compare(2,1) -> 10/1 
             meaning 2 is 2(twice) bigger than 1
        """
        if self.has_mathching_units(other):
            a = self.get_value()
            b = other.get_value()
            dif = abs(a/b)
            print(f"The size of {other} is {dif:.4f}% bigger than {self}")
            return dif
        else:
            raise Exception(f"Can't compare {self._u_text()} to {other._u_text()}")
    def __pow__(self, exponent):
        res = Unit(1,"m")
        res._value = self.get_value()**(exponent)
        res._units = [[i[0],i[1]*exponent] for i in self.get_units()]
        return res
    def __neg__(self):
        res = Unit(-1 * self._value, self._u_text())
        return res

    ### Deals with convertions of units & changes self
    def to(self, u_str):
        if u_str == "imperial":
            self.to_metric()
            u_list = [i[0] for i in self._units] # current units
            factor = 1
            for i,u in enumerate(u_list):
                self._units[i][0] = iconv[u] #Change the unit to imperial
                factor *= convertions[mconv[iconv[u]]][iconv[u]]**(-1 * self._units[i][1])
            self._value *= factor # Correct the value
            return self
        elif u_str == "metric":
            u_list = [i[0] for i in self._units] # current units
            factor = 1
            for i,u in enumerate(u_list):
                self._units[i][0] = mconv[u] #Change the unit to metric
                factor *= convertions[mconv[u]][u]**(self._units[i][1])
            self._value *= factor # Correct the value
            return self
        else:
            self.to_metric()
            u_list = u_str.split(",")
            s_list = [i[0] for i in self._units]
            factor = 1
            try: 
                for u in u_list:
                    index = s_list.index(mconv[u])
                    self._units[index][0] = u
                    factor *= convertions[mconv[u]][u]**(-1 * self._units[index][1])
                self._value *= factor
                return self
            except Exception as e:
                print(f"Some unit wasnt compatable with {self}\nFULL LOG:{e}")
                return None

class Constants(object):
    # From this sites:
    #   https://physics.info/constants/
    #   https://physics.nist.gov/cuu/pdf/all_2002.pdf
    pi = 3.1415926535879392                         # pi
    tau = 2*pi                                      # tau
    g = Unit(9.81, "m.s^-2")                        # Earth gravity
    c = Unit(299_792_458, "m.s^-1")                 # Speed of light in vacuum
    mu0 = Unit(4*pi*1e-7, "N.Amp^-2")               # Magnetic Constant
    e0 = Unit(8.8541878128e-12,"C^2.N^-1.m^-2")     # Electric/E Permitivity Constant
    K = e0**(-1) * (1/(4*pi))                       # Coulomb's Constant
    Z0 = Unit(376.730313461, "Ohm")                 # Vaccum impedance
    h = Unit(1.986445857e-25, "J.m")                # Planck's Constant
    hbar = h * (1/tau)                              # Hbar
    G = Unit(6.6742e-11, "m^3.kg^-1.s^-2")          # Gravitational Constant
    Q_electron = Unit(-1.602176e-19, "C")           # Charge of an electron
    m_eletron = Unit(9.1093826e-31, "kg")           # Mass of an electron
    m_proton = Unit(1.67262171e-27, "kg")           # Mass of a proton
    m_earth = Unit(5.9722e24, "kg")                 # Mass of the earth
    m_moon = Unit(7.348e22, "kg")                   # Mass of the moon
    m_sun = Unit(1.988500e30, "kg")                 # Mass of the sun
    k = Unit(1.380649e-23, "J.K^-1")                # Boltzmann's Constant
    R = Unit(8.314462618, "J.mol^-1.K^-1")          # Gas Constant
    sigma = Unit(5.670374419e-8, "W.m^-2.K^-4")     # Stefan-Boltz­mann's Constant
    H0 = Unit(2.25e-18, "s^-1")                     # Hubble's Constant
    a0 = Unit(5.29177210903e-11, "m")               # Bohr's radius
    alpha = 1/137                                   # Fine-structure constant

# TODO : implement this spiel
# import numpy as np
# class npUnit(Unit):
#     def __init__(self, value: np.ndarray, unit_str: str):


if __name__ == "__main__":
    u = Unit
    # hb = Constants.hbar
    # h = Constants.h
    # mu0 = Constants.mu0
    # test = u(1, "N^2.J^-2")
    # h.compare(hb)

    g = Constants.g
    velocity = u(1, "km.s^-1") 
    time = u(2, "min") 
    d1 = - g*time**2*0.5
    d2 = velocity*time

    distance = d1+ d2

    print( distance )                   # Unit( 2.88000e+04 m )  \n
    print( distance.to("au") )   # 28.80000 km 

    # print(Constants.K)
