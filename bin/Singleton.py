class Singleton(object):
# Adopted from http://stackoverflow.com/questions/6760685/creating-a-singleton-in-python
##  Metaclass implementation of a python singleton
#   __new__( ... )  : when a new Singleton object is called a new object will be created iff
#                     no instances of this object exists,
#                     otherwise a reference to that object is returned
    _instances = {}
    def __new__(class_, *args, **kwargs):
        if class_ not in class_._instances:
            class_._instances[class_] = super(Singleton, class_).__new__(class_, *args, **kwargs)
        return class_._instances[class_]
