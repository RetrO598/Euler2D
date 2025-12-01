#pragma once

#include <common/common.hpp>
#include <memory>
namespace preprocess {
class ElementConcept {
 public:
  virtual ~ElementConcept() = default;
  virtual std::unique_ptr<ElementConcept> clone() const = 0;
  virtual Index getnNodes() const = 0;
  virtual Index getnFaces() const = 0;
  virtual Index getnNeighbors() const = 0;
  virtual Index getNeighborElems(Index iface) const noexcept = 0;
  virtual void setNeighborElems(Index iface, Index elemIndex) noexcept = 0;
  virtual Index getNodes(Index inode) const noexcept = 0;
  virtual void setNodes(Index inode, Index nodeIndex) noexcept = 0;
  virtual Index getMaxNodesFace() const noexcept = 0;
  virtual Index getVTKtype() const noexcept = 0;
  virtual real getCoordCenter(Index idim) const noexcept = 0;
  virtual void setCoordCenter(Index idim, real icoord) noexcept = 0;
  virtual real getCoordFace(Index iface, Index idim) const noexcept = 0;
  virtual void setCoordFace(Index iface, Index idim, real icoord) noexcept = 0;
  virtual Index getNodesNeighbor(Index inode,
                                 Index ineighbor) const noexcept = 0;
  virtual Index getnNodesFace(Index iface) const noexcept = 0;
  virtual Index getnNeighborNodes(Index inode) const noexcept = 0;
  virtual Index getFaces(Index iface, Index inode) const noexcept = 0;
  virtual void changeOrientation() noexcept = 0;
};

template <class Type>
class ElementModel final : public ElementConcept {
 public:
  Type elem;

  ElementModel(const Type& type) : elem(type){};

  std::unique_ptr<ElementConcept> clone() const override {
    return std::make_unique<ElementModel<Type>>(*this);
  }
  Index getnNodes() const override { return Type::getnNodes(); }
  Index getnFaces() const override { return Type::getnFaces(); }
  Index getnNeighbors() const override { return Type::getnNeighbors(); }
  Index getNeighborElems(Index iface) const noexcept override {
    return elem.getNeighborElems(iface);
  }
  void setNeighborElems(Index iface, Index elemIndex) noexcept override {
    elem.setNeighborElems(iface, elemIndex);
  }
  Index getNodes(Index inode) const noexcept override {
    return elem.getNodes(inode);
  }
  void setNodes(Index inode, Index nodeIndex) noexcept override {
    elem.setNodes(inode, nodeIndex);
  }
  Index getMaxNodesFace() const noexcept override {
    return elem.getMaxNodesFace();
  }
  Index getVTKtype() const noexcept override { return elem.getVTKtype(); }
  real getCoordCenter(Index idim) const noexcept override {
    return elem.getCoordCenter(idim);
  }
  void setCoordCenter(Index idim, real icoord) noexcept override {
    elem.setCoordCenter(idim, icoord);
  }
  real getCoordFace(Index iface, Index idim) const noexcept override {
    return elem.getCoordFace(iface, idim);
  }
  void setCoordFace(Index iface, Index idim, real icoord) noexcept override {
    elem.setCoordFace(iface, idim, icoord);
  }
  Index getNodesNeighbor(Index inode, Index ineighbor) const noexcept override {
    return elem.getNodesNeighbor(inode, ineighbor);
  }
  Index getnNodesFace(Index iface) const noexcept override {
    return elem.getnNodesFace(iface);
  }
  Index getnNeighborNodes(Index inode) const noexcept override {
    return elem.getnNeighborNodes(inode);
  }
  Index getFaces(Index iface, Index inode) const noexcept override {
    return elem.getFaces(iface, inode);
  }
  void changeOrientation() noexcept override { elem.changeOrientation(); }
};

class ElementHandle {
 public:
  template <class T>
  ElementHandle(const T& elem)
      : self(std::make_unique<ElementModel<T>>(elem)) {}

  Index getnNodes() const { return self->getnNodes(); }
  Index getnFaces() const { return self->getnFaces(); }
  Index getnNeighbors() const { return self->getnNeighbors(); }
  Index getNeighborElems(Index iface) const noexcept {
    return self->getNeighborElems(iface);
  }
  void setNeighborElems(Index iface, Index elemIndex) noexcept {
    self->setNeighborElems(iface, elemIndex);
  }
  Index getNodes(Index inode) const noexcept { return self->getNodes(inode); }
  void setNodes(Index inode, Index nodeIndex) noexcept {
    self->setNodes(inode, nodeIndex);
  }
  Index getMaxNodesFace() const noexcept { return self->getMaxNodesFace(); }
  Index getVTKtype() const noexcept { return self->getVTKtype(); }
  real getCoordCenter(Index idim) const noexcept {
    return self->getCoordCenter(idim);
  }
  void setCoordCenter(Index idim, real icoord) noexcept {
    self->setCoordCenter(idim, icoord);
  }
  real getCoordFace(Index iface, Index idim) const noexcept {
    return self->getCoordFace(iface, idim);
  }
  void setCoordFace(Index iface, Index idim, real icoord) noexcept {
    self->setCoordFace(iface, idim, icoord);
  }
  Index getNodesNeighbor(Index inode, Index ineighbor) const noexcept {
    return self->getNodesNeighbor(inode, ineighbor);
  }
  Index getnNodesFace(Index iface) const noexcept {
    return self->getnNodesFace(iface);
  }
  Index getnNeighborNodes(Index inode) const noexcept {
    return self->getnNeighborNodes(inode);
  }
  Index getFaces(Index iface, Index inode) const noexcept {
    return self->getFaces(iface, inode);
  }
  void changeOrientation() noexcept { return self->changeOrientation(); }

 private:
  std::unique_ptr<ElementConcept> self;
};
}  // namespace preprocess