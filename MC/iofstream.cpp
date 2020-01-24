#include <iostream>

struct OutputAndConsole : std::ofstream
{
    OutputAndConsole(const std::string& fileName)
       : std::ofstream(fileName)
       , fileName(fileName)
    {};

    const std::string fileName;
};

template <typename T>
OutputAndConsole& operator<<(OutputAndConsole& os, const T& var)
{
	std::cout << (var);
	static_cast<std::ofstream&>(os)<< (var);
	return (os);
};

/*
template <typename T>
OutputAndConsole& operator<<(OutputAndConsole& strm, const T& var)
{
    std::cout << var;
    static_cast<std::ofstream&>(strm) << var;
    return strm;
};
*/
