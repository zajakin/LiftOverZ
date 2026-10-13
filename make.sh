  VER="1.0.0" 
  # git tag v${VER} && git push origin v${VER}
  # git tag -d v${VER} && git push origin --delete v${VER}
  git tag -a v${VER} -m "Release LiftOverZ v${VER}"
  git push origin main --tags
  