import re
import json


comparable = [
    'int',
    'uint',
    'unint_32t',
    'bool',
    'string',
    'std::string',
    'float',
    'short',
    'double',
    'long',
    'long long',
    'long double',
    'char *',
    'char*',
    'char',
    'Rect'
]

removable = [
    'unsigned',
    'static',
    'const'
]

def stem(s):
    src = s.split(' ')
    target = []
    for item in src:
        if item not in removable:
            target.append(item)
    return ' '.join([str(elem) for elem in target])

def getStruct(name, structs):
    for struct in structs:
        if(struct['name']==name):
            return struct
    return None

def components(structs, name, _type):
    if(stem(_type) in comparable or isRef(_type)):
        return [name]
    idx = 0;
    absType = _type.rstrip(' *')
    struct = getStruct(_type.rstrip(' *'),structs)
    if struct is not None:
        ret = []
        for field in struct['fields']:
            target = components(structs, field['name'], field['type'])
            if(_type.find('*')==-1):
                ret.extend([name+'.'+str(elem) for elem in target])
            else:
                ret.extend([name+'->'+str(elem) for elem in target])
        return ret
    return []

def addOnceToList(src, target):
    if isinstance(src, list): 
        for _obj in src:
            if _obj not in target:
                target.insert(0, _obj)
    elif src not in target:
        target.add(0, src)
    return target

def addOnceToDict(src, target):
    if isinstance(src, list): 
        for _obj in src:
            target.setdefault(_obj, [])
    elif src not in target:
        target.setdefault(src, [])
    return target

def isBitFields(field,structs):
    if 'bits' in field:
        return True
    struct = getStruct(field['type'], structs)
    if struct is None:
        return False
    for field in struct['fields']:
        if isBitFields(field, structs):
            return True
    return False
    
def getFunctionalName(name):
    if name.islower():
        return ''.join([n.capitalize() for n in name.split('_')])
    else:
        return name
    
def getClassIndex(schema, name):
    for i in range(len(schema['classes'])):
        if schema['classes'][i]['name']==name:
            return i
    return -1

def getTableName(name):
    if len(name)>2 and name[:2]=='db':
        name = name[2:]
    return '_{}_tbl'.format(name.lower())

def isRef(type_name):
    return True if type_name.startswith("dbId<")  and type_name[-1] == '>' else False

def isHashTable(type_name):
    return True if type_name.startswith("dbHashTable<")  and type_name[-1] == '>' else False

def getHashTableType(type_name):
    if not isHashTable(type_name) or len(type_name) < 13:
        return None

    return type_name[12:-1] + "*"

def isTemplateType(type_name):
    openBracket = type_name.find("<")
    if openBracket == -1:
        return False

    closedBracket = type_name.find(">")

    return False if closedBracket == -1 or closedBracket < openBracket else True

def getTemplateType(type_name):
    if isTemplateType(type_name) == False:
        return None 
    numBrackets = 1
    
    openBracket = type_name.find("<")
    for i in range(openBracket + 1, len(type_name)):
        
        if type_name[i] == "<":
            numBrackets += 1
        elif type_name[i] == ">":
            numBrackets -=1
            if numBrackets == 0:
                closedBracket = i
            
    return type_name[openBracket+1:closedBracket]



def getRefType(type_name):
    if not isRef(type_name) or len(type_name) < 7:
        return None

    return type_name[6:-1] + "*"
   
