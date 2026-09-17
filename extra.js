// Click a figure to see it large in front of the page; click again or press
// Esc to close. pkgdown loads this file on every page of the site.
// Only images shown smaller than their real size can be enlarged, so small
// icons are left alone, and images inside links keep their link.
document.addEventListener("DOMContentLoaded", function () {
  var overlay = document.createElement("div");
  overlay.className = "rvtr-lightbox";
  overlay.setAttribute("role", "dialog");
  overlay.setAttribute("aria-label", "Enlarged figure");
  var big = document.createElement("img");
  overlay.appendChild(big);
  document.body.appendChild(overlay);

  function close() {
    overlay.classList.remove("open");
    document.body.classList.remove("rvtr-lightbox-open");
  }
  overlay.addEventListener("click", close);
  document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") close();
  });

  function enlargeable(img) {
    return img.naturalWidth > img.clientWidth + 1;
  }
  function mark(img) {
    img.classList.toggle("rvtr-zoomable", enlargeable(img));
  }

  document.querySelectorAll("main img").forEach(function (img) {
    if (img.closest("a")) return;
    if (img.complete) mark(img); else img.addEventListener("load", function () { mark(img); });
    window.addEventListener("resize", function () { mark(img); });
    img.addEventListener("click", function () {
      if (!enlargeable(img)) return;
      big.src = img.currentSrc || img.src;
      big.alt = img.alt;
      overlay.classList.add("open");
      document.body.classList.add("rvtr-lightbox-open");
    });
  });
});
