#ifndef MC_EXCEPTIONS_H
#define MC_EXCEPTIONS_H

#include <exception>


#define EXIT(mcEx) { try{ throw mcEx; } catch(std::exception &e){ std::cerr << e.what() << std::endl; exit(EXIT_FAILURE); } }


class mcException : public std::exception {
    virtual const char *what() const throw() = 0;
};


class fellerCondition : public mcException {
    virtual const char *what() const throw() {
        return "Feller condition not respected.";
    }
};


#endif //MC_EXCEPTIONS_H
