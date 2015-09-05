#ifndef XYZ_GEOM_KDTREE_HPP_
#define XYZ_GEOM_KDTREE_HPP_

//
// K-Dimensional Tree (Axis Aligned Binary Space Division Tree)
//
namespace xyz {
class Intersection;
class Triangle;
class AbstractKDNode;

//
// kD-tree
//
class KDTree {
  HI_DISALLOW_COPY_AND_ASSIGN(KDTree);

 public:
  KDTree(__in std::size_t max_buket_size = 16,
         __in std::size_t max_hierarchy_depth = 16);
  ~KDTree();

 public:
  // property
  void SetMaxBucketSize(__in std::size_t max_bucket_size);
  void SetMaxHierarchyDepth(__in std::size_t max_hierarchy_depth);

 public:
  // builder and cleaner
  void Clear();
  void Build(__in std::vector<Triangle const *> const &set);

 private:
  AbstractKDNode *Build(__in std::vector<Triangle const *> const &set,
                        __in float3_t const &min, __in float3_t const &max,
                        __in std::size_t depth);

 public:
  // find Intersection
  Triangle const *FindIntersection(__in float3_t const &vOrigin,
                                   __in float3_t const &vDirection,
                                   __out Intersection *const pParam) const;

 private:
  std::size_t max_bucket_size_;
  std::size_t max_hierarchy_depth_;

  AbstractKDNode *root_;
  float3_t min_;
  float3_t max_;
};

}  // end of namespace xyz

#endif XYZ_GEOM_KDTREE_HPP_
