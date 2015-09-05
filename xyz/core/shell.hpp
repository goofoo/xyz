#ifndef XYZ_SHELL_HPP_
#define XYZ_SHELL_HPP_

namespace xyz {
class Scene;

std::tstring const &GetOutdirPath();

int GetRenderingTime();
int GetIntervalTime();

bool IsFinished();
void FinishRendering();

void Display();
void Resize(int, int);
void Keyboard(unsigned char, int, int);

int Shell(std::tistream &in, bool const gui = true);

}  // end of namespace xyz

#endif XYZ_SHELL_HPP_
