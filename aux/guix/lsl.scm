;; TODO: There is not much missing for this to be submitted as a patch to Guix
;; See `https://guix.gnu.org/manual/en/html_node/Submitting-Patches.html`

;;(define-module (gnu packages lsl)
(define-module (lsl)
  #:use-module (guix packages)
  #:use-module (guix git-download)
  #:use-module (guix build-system cmake)
  #:use-module ((guix licenses)
                #:prefix license:)
  #:use-module (gnu packages xml))

(define-public lsl
  (package
    (name "lsl")
    (version "v1.16.2")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/sccn/liblsl")
             (commit version)))
       (file-name (git-file-name name version))
       (sha256
        (base32 "1six8qd65m3q03wbi2fyd3i1bpl9g2ag97f37c39nlrq34mvnswy"))))
    (build-system cmake-build-system)
    (arguments
     '(#:configure-flags (list "-DLSL_UNIXFOLDERS=ON"
                               "-DLSL_BUNDLED_PUGIXML=OFF")
       #:tests? #f))
    (inputs (list pugixml))
    (home-page "https://labstreaminglayer.org")
    (synopsis "Lab Streaming Layer")
    (description
     "LSL is an open-source networked middleware ecosystem to
stream, receive, synchronize, and record neural, physiological, and behavioral
data streams acquired from diverse sensor hardware.")
    (license license:expat)))

lsl
