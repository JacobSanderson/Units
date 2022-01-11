# Another related proyect
#  - https://scimath.readthedocs.io/en/latest/index.html
import numpy as np

metricu = "g,m,s,C,mol,K".split(",")
metric_order = {i:j for i,j in  zip( metricu, range(len(metricu)) ) }

# the prefixes in units are
prefixes = {
 "P"  : 1e15,    # peta
 "T"  : 1e12,    # tera
 "G"  : 1e9,     # giga
 "M"  : 1e6,     # mega
 "k"  : 1e3,     # kilo
 "d"  : 1e-2,    # deci
 "m"  : 1e-3,    # mili
 "u"  : 1e-6,    # micro
 "n"  : 1e-9,    # nano
 "p"  : 1e-12,   # pico
 "f"  : 1e-15,   # femto
}

# equivalences to SI/metric, e.g.  { 'mi' : 1.6e3 (meters) , ... }
convertions = {
    "s"  : {                
        "s"  : 1,           # Second           
        "year":3.15576e7,   # Year
        "day": 86400,       # Day
        "hr" : 3600,        # Hour
        "min": 60,          # Minute
    },
    "g" : {      
        "g"  : 1,           # gram
        "ton": 9.071847e2,  # ton
        "lb" : 0.4535924,   # pound
    },
    "m"  : {     
        "m"  : 1,           # Meter
        "A"  : 1e-10,       # Angstrom
        "in" : 0.0254,      # Inches
        "ft" : 0.3048,      # Feet
        "yd" : 0.9144,      # Yard
        "mi" : 1.609344e3,  # Mile
        "au" : 1.495979e11  # Astronomical unit
    },
    "mol" : {
        "mol" : 1,          # Mol
        "particles" : 1/(6.02214076e23) # Particles
    },
    "C" : {    
        "C" : 1             #Coulomb
    },
    "K" : {             
        "K" : 1             #Kelvin             
    }
}

all_units = []
for i in convertions.keys():
    for j in convertions[i].keys():
        all_units.append(j)

# cu := Combined Unit
# fu := Fundamental Unit
#   cu   : (n: "1 cu = n fu" , fu's )
combined_units = {
    "V"  : (1,"kg.m^2.s^-2.C^-1"),    # Volts
    "Amp": (1,"s^-1.C"),              # Ampere
    "Ohm": (1,"kg.m^2.s^-1.C^-2"),    # Ohm
    "N"  : (1,"kg.m.s^-2"),           # Newton
    "J"  : (1,"kg.m^2.s^-2"),         # Joule
    "cal": (4.184,"kg.m^2.s^-2"),     # Calorie
    "Pa" : (1,"kg.m^-1.s^-2"),        # Pascal
    "PSI": (6894.757,"kg.m^-1.s^-2"), # Pounds per Square Inch
    "ba" : (100_000,"kg.m^-1.s^-2"),  # Bars
    "atm": (101_325,"kg.m^-1.s^-2"),  # atmfosfere
    "L"  : (1,"m^3"),                 # Liters
    "ga" : (3.785,"m^3"),             # Gallons
    "W"  : (1,"kg.m^2.s^-3")          # Watts
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
        self._value = value
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
                    cu_list = self.get_u_matrix(combined_units[cu][1])
                    cu_list = [ [i[0],exponent*i[1]] for i in cu_list]
                    for i in cu_list:
                        self._value *= (combined_units[cu][0])**(exponent*i[1])
                    temp = temp.replace(cu+f"^{exponent}", self._u_text2(cu_list))
                else:
                    temp = temp.replace(cu,combined_units[cu][1])
                    self._value *= combined_units[cu][0] # Warning : Not sure if is correct
        self._units = self.get_u_matrix(temp)
        # self._units has the form [["kg",1],["hr",4],...]
        # may be repeated, cuz of <<Combined Units>> so
        self._s_simplify()      # Eliminate repeated values
        self.to_metric()        # Convert to metric
        self._units.sort(key = lambda s:metric_order[s[0]]) # Just for kiks

    ### Utility functions
    def get_u_matrix(self, u_str:str) -> list:
        """Input:  "kg^2.s^-3.cm^5" | Ouput: [['kg',2],['s',-3],['cm',5]]"""
        t_input = u_str.split(".")
        res = []
        for i,u in enumerate(t_input):
            if "^" in u:
                index = u.rfind("^")
                exponent = int( u[index+1:] )
                unit = u[:index]
            else:
                unit = u
                exponent = 1
            res.append( [unit, exponent] )
            try:
                # Checking for prefixes
                if unit[0] in prefixes.keys():
                    # Avoid mixing prefixes and units up
                    if unit[0] in "mdfup":  
                        isdiff = True
                        for c_unit in all_units:
                            if unit == c_unit:
                                isdiff = False
                                res[i] = [unit, exponent]
                                break
                        if isdiff:
                            res[i] = [unit[1:], exponent]
                            self._value *= prefixes[unit[0]]
                    else:
                        res[i] = [unit[1:], exponent]
                        self._value *= prefixes[unit[0]]
                else:
                    res[i] = [unit, exponent]
            except Exception as e:
                print(f"ERROR : Couldn't parse {u_str}, maybe is not written correctly.\n{e}")
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
        return "Unit( {:.5e} {:s} )".format(self.get_value(), self._u_text())
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
        if type(other) == float:
            return Unit(self.get_value()/other, self._u_text())
        else:
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
            dummy = "bigger"
            if dif < 1:
                dummy = "smaller"
                dif = 1/dif
            print(f"The size of {other} is {dif:.4} times {dummy} than {self}")
            return dif
        else:
            raise Exception(f"Can't compare {self._u_text()} to {other._u_text()}")
    def __pow__(self, exponent):
        res = Unit(1,"m")
        res._value = self.get_value()**(exponent)
        res._units = [[i[0],i[1]*exponent] for i in self.get_units()]
        return res
    def __neg__(self):
        return Unit(-1.0 * self._value, self._u_text())

    ### Deals with convertions of units & CHANGES self !!!!
    def to(self, u_str:str):
        try:
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
                u_list = u_str.split(".")            #["hr","ft"]
                s_list = [i[0] for i in self._units] #["m","s","C"]
                factor = 1
                for u in u_list:
                    index = s_list.index(mconv[u])   #mconv: unit -> metric unit
                    self._units[index][0] = u        #Change the name of units
                    factor *= convertions[mconv[u]][u]**(-1 * self._units[index][1])
                self._value *= factor 
                return self
        except Exception as e:
            print("\n")
            print(f"Error getting from {self.get_units()} to {u_str}.\nFULL LOG:{e}")

# From this sites:
#   https://physics.info/constants/
#   https://physics.nist.gov/cuu/pdf/all_2002.pdf
class Constant(Unit):
    def __init__(self, value: float, unit_str: str, name: str):
        super().__init__(value, unit_str)
        self._name = name
    def __repr__(self) -> str:
        return f"{self._name} : {self._value} {self._u_text()} "

pi = 3.1415926535879392                                                # pi
tau = 2*pi                                                             # tau
g = Constant(9.81, "m.s^-2", "Earth gravity")                          # Earth gravity
c = Constant(299_792_458, "m.s^-1", "Speed of light in vacuum")        # Speed of light in vacuum
mu0 = Constant(4*pi*1e-7, "N.Amp^-2", "Magnetic Constant")             # Magnetic Constant
e0 = Constant(8.8541878128e-12,"C^2.N^-1.m^-2", "Electric Constant")   # Electric/E Permitivity Constant
K = Constant((e0**(-1)*(1/(4*pi))).n(),"C^2.N^-1.m^-2","Coulomb's Constant")#  Coulomb's Constant
Z0 = Constant(376.730313461, "Ohm", "Vaccum impedance")                # Vaccum impedance
h = Constant(1.986445857e-25, "J.m", "Planck's Constant")              # Planck's Constant
hbar = Constant( (h/tau).n() , "J.m", "Hbar")                          # Hbar
G = Constant(6.6742e-11, "m^3.kg^-1.s^-2", "Gravitational Constant")   # Gravitational Constant
Q_electron = Constant(-1.602176e-19, "C", "Charge of an electron")     # Charge of an electron
m_electron = Constant(9.1093826e-31, "kg", "Mass of an electron")      # Mass of an electron
m_proton = Constant(1.67262171e-27, "kg", "Mass of a proton")          # Mass of a proton
m_earth = Constant(5.9722e24, "kg", "Mass of the earth")               # Mass of the earth
m_moon = Constant(7.348e22, "kg", "Mass of the moon")                  # Mass of the moon
m_sun = Constant(1.988500e30, "kg", "Mass of the sun")                 # Mass of the sun
k = Constant(1.380649e-23, "J.K^-1", "Boltzmann's Constant")           # Boltzmann's Constant
R = Constant(8.314462618, "J.mol^-1.K^-1", "Gas Constant")             # Gas Constant
sigma = Constant(5.670374419e-8, "W.m^-2.K^-4", "Stefan-Boltz­mann's Constant")     # Stefan-Boltz­mann's Constant
H0 = Constant(2.25e-18, "s^-1", "Hubble's Constant")                   # Hubble's Constant
a0 = Constant(5.29177210903e-11, "m", "Bohr's radius")                 # Bohr's radius
alpha = 1/137                                                          # Fine-structure constant
__map_dic = {
    "pi" : pi,
    "tau" : tau,
    "g" : g,
    "earth gravity" : g,
    "c" : c,
    "speed light" : c,
    "speed of light" : c,
    "mu0" : mu0,
    "magnetic" : mu0,
    "e0" : e0,
    "electric permitivity" : e0,
    "coulomb" : K,
    "z0" : Z0,
    "vaccum impedance" : Z0,
    "h" : h,
    "planck" : h,
    "hbar" : hbar,
    "gravitational" : G,
    "charge electron" : Q_electron,
    "charge of an electron" : Q_electron,
    "mass electron" : m_electron,
    "mass of an electron" : m_electron,
    "mass proton" : m_proton,
    "mass of a proton" : m_proton,
    "mass earth" : m_earth,
    "mass of the earth" : m_earth,
    "mass moon" : m_moon,
    "mass of the moon" : m_moon,
    "mass sun" : m_sun,
    "mass of the sun" : m_sun,
    "boltzmann" : k,
    "gas" : R,
    "ideal gas" : R,
    "stefan boltz­mann" : sigma,
    "stefan-boltz­mann" : sigma,
    "hubble" : H0,
    "bohr radius" : a0,
    "bohr's radius" : a0,
    "fine structure" : alpha,
    "fine-structure" : alpha,
}
def get_constant(name:str) -> Constant:
    try:
        if name == "G":
            return G
        else:
            return __map_dic[name.lower()]
    except Exception as e:
        print(f"There is no constant named {name} \n",e)
    

# TODO : implement this spiel
# import numpy as np
# class npUnit(Unit):
#     def __init__(self, value: np.ndarray, unit_str: str):

if __name__ == "__main__":
    def test_00():
        # print(K)
        a = h.compare(hbar)
        print(a)
    def test_02():
        a = get_constant("speed of light")
        b = get_constant("speed light")
        c = get_constant("c").to("km")
        print(a,b,c, sep="\n")

    def test_01():
        g = get_constant("g")
        velocity = Unit(1, "km.s^-1")
        time = Unit(2, "min")
        distance = - g*time**2*(1/2) + velocity * time * 2

        print( distance )                   # Unit( 2.88000e+04 m )  \n
        print( distance.to("mi").text() )   # 28.80000 km
    
    def test_03():
        a = Unit(1, "V")
        b = Unit(1, "Amp")
        c = Unit(1, "s")
        res = a*b*c

        # Quiero que esto funcione:
        # print( res.to("kV") )
    test_00()
