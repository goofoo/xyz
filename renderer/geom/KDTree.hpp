#ifndef __TGIR_GEOM_KDTREE_HPP__
#define __TGIR_GEOM_KDTREE_HPP__

//
// K-Dimensional Tree (Axis Aligned Binary Space Division Tree)
//
namespace tgir {
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
  void Build(__in std::vector<tgir::Triangle const *> const &set);

 private:
  AbstractKDNode *Build(__in std::vector<tgir::Triangle const *> const &set,
                        __in tgir::Vector3 const &min,
                        __in tgir::Vector3 const &max, __in std::size_t depth);

 public:
  // find Intersection
  tgir::Triangle const *FindIntersection(
      __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
      __out tgir::Intersection *const pParam) const;

 private:
  std::size_t max_bucket_size_;
  std::size_t max_hierarchy_depth_;

  tgir::AbstractKDNode *root_;
  tgir::Vector3 min_;
  tgir::Vector3 max_;
};

}  // end of namespace tgir

#endif __TGIR_GEOM_KDTREE_HPP__
