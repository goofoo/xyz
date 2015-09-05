// main.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//
#include "stdafx.hpp"

static long nvertices;
static long ntriangles;
static long log_nvertices;
static long log_ntriangles;
static long ivertices;
static long itriangles;

std::vector<int> tristrips;

static int vertex_cb(p_ply_argument argument) {
  long index;
  FILE *fp;
  ply_get_argument_user_data(argument, (void **)&fp, &index);

  switch (index) {
  case 0:
    if (ivertices == 0) {
      fprintf(fp, "# vertex %ld\n", nvertices);
      fprintf(stderr, "  vertices (%ld):   0%\r", nvertices);
    }
    fprintf(fp, "v %g", -ply_get_argument_value(argument));  // x軸だけ反転
    break;
  case 1:
    fprintf(fp, " %g ", ply_get_argument_value(argument));
    break;
  case 2:
    fprintf(fp, "%g\n", ply_get_argument_value(argument));
    if (++ivertices % log_nvertices == 0) {
      fprintf(stderr, "  vertices (%ld): % 3ld%%\r", nvertices,
              (100 * ivertices) / nvertices);
    }
    if (ivertices >= nvertices) {
      fputs(
          ""
          "\n",
          fp);
      fprintf(stderr, "  vertices (%ld): 100%\n", nvertices);
    }
    break;
  }

  return 1;
}

static int face_cb(p_ply_argument argument) {
  FILE *fp;
  ply_get_argument_user_data(argument, (void **)&fp, NULL);

  long length;
  long value_index;

  ply_get_argument_property(argument, NULL, &length, &value_index);

  if (value_index < 0) {
    if (itriangles == 0) {
      fprintf(fp, "# face %ld\n", ntriangles);
      fprintf(stderr, "  triangles (%ld):   0%\r", ntriangles);
    }
    fprintf(fp, "f %d", length);
  } else if (value_index == (length - 1)) {
    fprintf(fp, "%g\n", ply_get_argument_value(argument) + 1);

    ++itriangles;

    if (itriangles % log_nvertices == 0) {
      fprintf(stderr, "  triangles (%ld): % 3ld%%\r", ntriangles,
              (100 * itriangles) / ntriangles);
    }

    if (itriangles >= ntriangles) {
      fputs(
          ""
          "\n",
          fp);
      fprintf(stderr, "  triangles (%ld): 100%\n", ntriangles);
    }
  } else {
    fprintf(fp, "%g ", ply_get_argument_value(argument) + 1);
  }

  return 1;
}

static int tristrips_cb(p_ply_argument argument) {
  long length;
  long value_index;

  ply_get_argument_property(argument, NULL, &length, &value_index);

  if (value_index < 0) {
    tristrips.clear();
    tristrips.resize(length);
  } else {
    tristrips[value_index] = (int)ply_get_argument_value(argument) + 1;

    if ((value_index + 1U) == tristrips.size()) {
      // 面数を計算して出力
      FILE *fp;
      ply_get_argument_user_data(argument, (void **)&fp, NULL);

      long nfaces = 0;
      // 数えるだけ
      {
        int nvrt = 0;
        for (std::size_t i = 0, size = tristrips.size(); i < size; ++i) {
          if (tristrips[i] <= 0) {
            // reset
            nvrt = 0;
          } else if (nvrt < 2) {
            // skip
            nvrt++;
          } else {
            nfaces++;
          }
        }
      }
      // 実際に出力
      fprintf(fp, "# face %ld\n", nfaces);
      {
        int nvrt = 0;
        bool flip = false;
        for (std::size_t i = 0, size = tristrips.size(); i < size; ++i) {
          if (tristrips[i] <= 0) {
            // reset
            nvrt = 0;
            flip = false;
          } else if (nvrt < 2) {
            // skip
            nvrt++;
          } else {
            if (flip) {
              fprintf(fp, "f %d %d %d\n", tristrips[i - 2], tristrips[i - 1],
                      tristrips[i]);
            } else {
              fprintf(fp, "f %d %d %d\n", tristrips[i - 1], tristrips[i - 2],
                      tristrips[i]);
            }
            flip = !flip;
          }
        }
      }
      fputs(
          ""
          "\n",
          fp);
    }
  }

  return 1;
}

static int convert(char filename[], FILE *fp) {
  p_ply ply = ply_open(filename, NULL);
  if (!ply) {
    fprintf(stderr, "convert error: can not open file `%s\'\n", filename);
    return 1;
  }

  if (!ply_read_header(ply)) {
    fprintf(stderr, "convert error: can not read header `%s\'\n", filename);
    return 1;
  }

  nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, fp, 0);
  ply_set_read_cb(ply, "vertex", "y", vertex_cb, fp, 1);
  ply_set_read_cb(ply, "vertex", "z", vertex_cb, fp, 2);

  ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", face_cb, fp, 0);
  ply_set_read_cb(ply, "tristrips", "vertex_indices", tristrips_cb, fp, 0);

  log_nvertices = nvertices / 100;
  if (log_nvertices < 1) log_nvertices = 1;
  log_ntriangles = ntriangles / 100;
  if (log_ntriangles < 1) log_ntriangles = 1;
  ivertices = 0;
  itriangles = 0;

  fputs(
      "g obj1"
      "\n",
      fp);

  if (!ply_read(ply)) {
    fprintf(stderr, "convert error: can not read data `%s\'\n", filename);
    return 1;
  }

  ply_close(ply);

  return 0;
}

int main(int argc, char *argv[]) {
  char filename[512];

  for (int i = 1; i < argc; ++i) {
    sprintf(filename, "%s.obj", argv[i]);

    FILE *fp = fopen(filename, "wt");
    if (fp) {
      fprintf(stderr, "convert `%s\' to `%s\'\n", argv[i], filename);
      convert(argv[i], fp);
      fclose(fp);
    }
  }

  return 0;
}
