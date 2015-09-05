#pragma once

#include <hi/lang.hpp>

namespace hi {
// DIBSection���g�������ڑ���ł���摜
class memory_image sealed {
 public:
  memory_image();

  bool create(hi::sint width, hi::sint height);
  void destroy();

  void clear();
  bool load(std::tchar_t const *filename);
  bool copy(memory_image const &image);
  void raster(HDC hdc) const;

  void swap(memory_image &image);

  inline bool valid() const;
  inline hi::sint width() const;
  inline hi::sint height() const;
  inline hi::dword &operator[](std::size_t i);
  inline hi::dword operator[](std::size_t i) const;
  inline hi::dword &operator()(std::size_t x, std::size_t y);
  inline hi::dword operator()(std::size_t x, std::size_t y) const;

 private:
  bool init(hi::sint width, hi::sint height);

 private:
  hi::sint width_;
  hi::sint height_;
  hi::dword *pixels_;

  // �n���h���֌W(�V�X�e���Ɉˑ����镔��)
  HDC hDc_;  ///< �r�b�g�}�b�v�֏������ނ��߂̃f�o�C�X�R���e�L�X�g�ւ̃n���h��
  HBITMAP hDib_;     ///< �r�b�g�}�b�v�ւ̃n���h��
  HGDIOBJ hOldDib_;  ///< �f�o�C�X�R���e�L�X�g���O�Ɏ����Ă����r�b�g�}�b�v
};

inline bool memory_image::valid() const { return pixels_ != nullptr; }

inline hi::sint memory_image::width() const { return width_; }

inline hi::sint memory_image::height() const { return height_; }

inline hi::dword &memory_image::operator[](std::size_t i) { return pixels_[i]; }

inline hi::dword memory_image::operator[](std::size_t i) const {
  return pixels_[i];
}

inline hi::dword &memory_image::operator()(std::size_t x, std::size_t y) {
  return pixels_[y * width() + x];
}

inline hi::dword memory_image::operator()(std::size_t x, std::size_t y) const {
  return pixels_[y * width() + x];
}

}  // end of namespace hi
