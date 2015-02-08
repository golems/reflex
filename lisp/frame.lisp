;;;; -*- mode: lisp -*-
;;;;
;;;; Copyright (c) 2013, Georgia Tech Research Corporation
;;;; Copyright (c) 2015, Rice University
;;;; All rights reserved.
;;;;
;;;; Author(s): Neil T. Dantam <ntd@gatech.edu>
;;;; Georgia Tech Humanoid Robotics Lab
;;;; Under Direction of Prof. Mike Stilman <mstilman@cc.gatech.edu>
;;;;
;;;;
;;;; This file is provided under the following "BSD-style" License:
;;;;
;;;;
;;;;   Redistribution and use in source and binary forms, with or
;;;;   without modification, are permitted provided that the following
;;;;   conditions are met:
;;;;
;;;;   * Redistributions of source code must retain the above copyright
;;;;     notice, this list of conditions and the following disclaimer.
;;;;
;;;;   * Redistributions in binary form must reproduce the above
;;;;     copyright notice, this list of conditions and the following
;;;;     disclaimer in the documentation and/or other materials provided
;;;;     with the distribution.
;;;;
;;;;   * Neither the name of copyright holder the names of its
;;;;     contributors may be used to endorse or promote products
;;;;     derived from this software without specific prior written
;;;;     permission.
;;;;
;;;;   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
;;;;   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
;;;;   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;;;   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
;;;;   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
;;;;   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;;;   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;;;   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
;;;;   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
;;;;   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
;;;;   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
;;;;   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
;;;;   POSSIBILITY OF SUCH DAMAGE.

(in-package :reflex)

;;; TODO: separate varying and fixed frames

(defmacro def-frame-var (var doc)
  (let ((frame (gensym "FRAME")))
    `(defun ,(intern (concatenate 'string "FRAME-" (string var)))
         (,frame)
       ,doc
       (destructuring-bind (,frame &key ,var &allow-other-keys) ,frame
         (assert (eq :frame ,frame))
         ,var))))

(def-frame-var name "name of the frame")
(def-frame-var type "type of the frame, fixed or revolute")
(def-frame-var parent "parent of this frame")
(def-frame-var axis "frame axis of rotation, (x y z)")
(def-frame-var configuration "name of frame configuration")
(def-frame-var quaternion "fixed transform quaternion")
(def-frame-var translation "translation (x y z)")
(def-frame-var offset "Configuration offset")

(defun make-frame (name type parent &key
                                      quaternion
                                      translation
                                      configuration
                                      axis
                                      offset)
  (let ((name (string-upcase name))
        (type (etypecase type
                (string
                 (let ((s (string-downcase type)))
                   (cond
                     ((string= "revolute" s)
                      :revolute)
                     ((string= "prismatic" s)
                      :prismatic)
                     ((string= "fixed" s)
                      :fixed)
                     (t
                      (error "Unknown joint type ~A" s)))))
                (symbol type))))
    `(:frame :name ,name
             :type ,type
             :parent ,(when parent (string-upcase parent))
             :axis ,axis
             :configuration ,(cond (configuration configuration)
                                   ((eq type :fixed) nil)
                                   (t (concatenate 'string "Q_" name)))
             :quaternion ,quaternion
             :translation ,translation
             :offset ,offset)))



(defparameter *test-frames*
  `((frame :name frame_0
           :configuration q0)
    (frame :name frame_1
           :configuration q1)
    (frame :name frame_2
           :configuration q2)))

(defun prefix-frame-1 (prefix frame parent)
  (labels ((prefix (val) (when val (concatenate 'string prefix (string val)))))
    (append (list :frame
                  :parent (if (and (frame-parent frame)
                                   (not (equal parent (frame-parent frame))))
                              (prefix (frame-parent frame))
                              parent))
            (loop for k = (cdr frame) then (cddr k)
               while k
               for key = (first k)
               for val = (second k)
               unless (eq key :parent)
               append (list key
                            (cond
                              ((or (eq :name key) (eq :configuration key))
                               (prefix val))
                              (t
                               val)))))))

(defun prefix-frames (parent prefix frames)
  (loop for f in frames
       collect (prefix-frame-1 prefix f parent)))


(defun make-fixed-frame (name parent orientation
                         &key
                         (x 0)
                         (y 0)
                         (z 0))
  `(:frame :name ,name
          :parent ,parent
          :type :fixed
          :translation (,x ,y ,z)
          :quaternion ,(let ((q (amino:quaternion orientation)))
                            (list (amino::quaternion-x q)
                                  (amino::quaternion-y q)
                                  (amino::quaternion-z q)
                                  (amino::quaternion-w q)))))

(defun emit-indices (stream indices max)
  (format stream "~{~&#define ~{~A ~D~}~}"
          (if max
              (cons (list max (length indices)) indices)
              indices)))

(defun emit-name-indices (frames &key (stream t) max)
  (emit-indices stream
                (loop
                   for frame in frames
                   for i from 0
                   collect (list (frame-name frame)
                                 i))
                max))


(defun emit-config-indices (frames &key (stream t) max)
  (let ((vars (remove-duplicates (loop
                                    for frame in frames
                                    for config = (frame-configuration frame)
                                    when config
                                    collect config)
                                 :test #'equal)))
    (emit-indices stream
                  (loop
                     for i from 0
                     for v in vars
                     collect (list v i))
                  max)))

(defun float-string (f)
  (if (numberp f)
      (map 'string (lambda (c)
                     (if (or (eq c #\d)
                             (eq c #\D))
                         #\E
                         c))
           (format nil "~F" f))
      f))

(defun emit-fixed-frame (frame stream e)
  (let ((q (frame-quaternion frame))
        (v (frame-translation frame))
        (ptr (format nil "(~A+ldE*~A)" e (frame-name frame))))
    (dotimes (i 4)
      (format stream "~&    ~A[AA_TF_QUTR_Q + ~D] = ~A;"
              ptr i (float-string (elt q i))))
    (dotimes (i 3)
      (format stream "~&    ~A[AA_TF_QUTR_T + ~D] = ~A;"
              ptr i (float-string (elt v i))))))

(defun emit-quat (stream ptr angle x y z)
  (let* ((n (sqrt (+ (* x x) (* y y) (* z z))))
         (x (/ x n))
         (y (/ y n))
         (z (/ z n)))
    (cond
      ((= x 1)  (format stream "~&    aa_tf_xangle2quat( ~A, ~A );" angle ptr))
      ((= x -1) (format stream "~&    aa_tf_xangle2quat( -~A, ~A );" angle ptr))
      ((= y 1)  (format stream "~&    aa_tf_yangle2quat( ~A, ~A );" angle ptr))
      ((= y -1) (format stream "~&    aa_tf_yangle2quat( -~A, ~A );" angle ptr))
      ((= z 1)  (format stream "~&    aa_tf_zangle2quat( ~A, ~A );" angle ptr))
      ((= z -1) (format stream "~&    aa_tf_zangle2quat( -~A, ~A );" angle ptr))
       (t
        (format stream "~&    {")
        (format stream "~&        static const double axis[3] = {~A,~A,~A};"
                (float-string x)
             (float-string y)
             (float-string z))
     (format stream "~&        aa_tf_axang2quat2( axis, ~A, ~A );"
             angle
             ptr)
     (format stream "~&    }")))))

(defun offset-config (phi offset)
  (if offset
      (format nil "(~A+~A)" (float-string phi) (float-string offset))
      (float-string phi)))

(defun emit-revolute-frame (frame stream q-array e-array)
  (assert (frame-configuration frame) ()
          "No configuration variable for frame ~A" (frame-name frame))
  (let ((v (frame-translation frame))
        (axis (frame-axis frame))
        (ptr (format nil "(~A + ldE*~A + AA_TF_QUTR_Q)" e-array (frame-name frame)))
        (phi (format nil "~A[~A*incQ]" q-array (frame-configuration frame)))
        (offset (frame-offset frame)))
    (emit-quat stream ptr (offset-config phi offset)
               (elt axis 0) (elt axis 1) (elt axis 2))
    (dotimes (i 3)
      (format stream "~&    ~A[AA_TF_QUTR_T + ~D] = ~A;"
              ptr i (float-string (elt v i))))))


(defun emit-prismatic-frame (frame stream q-array e-array)
  (assert (frame-configuration frame) ()
          "No configuration variable for frame ~A" (frame-name frame))
  (let ((v (frame-translation frame))
        (axis (frame-axis frame))
        (ptr (format nil "(~A + ldE*~A + AA_TF_QUTR_Q)" e-array (frame-name frame)))
        (q (frame-quaternion frame))
        (phi (format nil "~A[~A*incQ]" q-array (frame-configuration frame)))
        (offset (frame-offset frame)))
    (dotimes (i 4)
      (format stream "~&    ~A[AA_TF_QUTR_Q + ~D] = ~A;"
              ptr i (float-string (elt q i))))
    (dotimes (i 3)
      (format stream "~&    ~A[AA_TF_QUTR_T + ~D] = ~A~A;"
              ptr i
              (float-string (elt v i))
              (cond ;;
                ((zerop (elt axis i)) "")
                ((= 1 (elt axis i))
                 (format nil " + ~A"
                         (offset-config phi offset)))
                (t
                 (format nil " + ~A * ~A"
                         (offset-config phi offset)
                         (float-string (elt axis i)))))))))

(defun emit-parents-array (name frames &key (stream t))
  (format stream "~&const ssize_t ~A[] = {~&~{~&    ~A~^,~}~&};"
          name
          (loop for f in frames
             for p = (frame-parent f)
             collect (if p p -1))))

(defun emit-axes-array (name frames &key (stream t))
  (format stream "~&const double ~A[][3] = {" name )
  (loop for f in frames
     for a = (map 'list #'float-string (frame-axis f))
     when a
     do (format stream "~&    {~{~A~^, ~}}, /* ~A */"
                a (frame-name f)))
  (format stream "~&};"))


(defun emit-names-array (name frames &key (stream t))
  (format stream "~&const char * ~A[] = {~&~{~&    \"~A\"~^,~}~&};"
          name
          (loop for f in frames
             for p = (frame-name f)
             collect (if p p -1))))


(defun rel-fun-decl (name &key (stream t) block-arrays)
  (format stream "~&void ~A( ~A )"
          name
          (if block-arrays
              "const double *AA_RESTRICT q, size_t incQ, double *AA_RESTRICT e, size_t ldE"
              "const double *AA_RESTRICT q, double *AA_RESTRICT e")))

(defun emit-rel-fun (name frames &key (stream t) block-arrays)
  (rel-fun-decl name :stream stream :block-arrays block-arrays)
  (format stream "~&{~&")
  (unless block-arrays
    (format stream "~&    const size_t incQ = 1;")
    (format stream "~&    const size_t ldE = 7;"))
  (loop
     for frame in frames
     for i from 0
     do
       (format stream "~&    // Frame: ~A, type: ~A"
               (frame-name frame)
               (frame-type frame))
       (case (frame-type frame)
         (:revolute (emit-revolute-frame frame stream "q" "e"))
         (:prismatic (emit-prismatic-frame frame stream "q" "e"))
         (:fixed (emit-fixed-frame frame stream "e"))
         (otherwise (error "Unknown type of frame ~A" (frame-name frame)))))
  (format stream "~&}"))

(defun frame-fun-decl (name &key (stream t))
  (format stream "~&void ~A( ~A )"
          name
          (concatenate 'string
                       "const double *AA_RESTRICT q, size_t incQ, "
                       "double *AA_RESTRICT E_rel, size_t ldRel, "
                       "double * AA_RESTRICT E_abs, size_t ldAbs,"
                       "int options")))

(defun emit-frame-fun (name rel-fun abs-fun frames &key (stream t))
  (frame-fun-decl name :stream stream)
  ;; TODO: interleave relative and absolute frame computation for more linear memory access
  ;; TODO: add options to elide recopying fixed frames
  (format stream "~&{~&")
  (format stream "~&    (void) options;")
  (format stream "~&    if(q) ~A(q, incQ, E_rel, ldRel);" rel-fun)
  (format stream "~&    if(E_abs) ~A(E_rel, ldRel, E_abs, ldAbs);" abs-fun)
  (format stream "~&}"))

;; (defun emit-jacobian-fun (name frames &key (stream t))
;;   (format stream "~&void ~A~&(~{~A~^, ~}) ~&{~&"
;;           name
;;           '("const double *restrict tf_abs"
;;             "size_t n" "const size_t *restrict indices" "const double *restrict pe"
;;             "double *restrict J" "size_t ldJ"))

(defun abs-fun-decl (name &key (stream t) block-arrays)
  (format stream "~&void ~A( ~A )"
          name
          (if block-arrays
              "const double * AA_RESTRICT rel, size_t ldRel, double * AA_RESTRICT abs, size_t ldAbs"
              "const double * AA_RESTRICT rel, double * AA_RESTRICT abs")))


(defun emit-abs-fun (name frames &key normalize (stream t) block-arrays)
  (abs-fun-decl name :stream stream :block-arrays block-arrays)
  (format stream "~&{~&")
  (unless block-arrays
    (format stream "~&    const size_t ldRel = 7;")
    (format stream "~&    const size_t ldAbs = 7;"))
  (loop
     for f in frames
     for name = (frame-name f)
     for rel = (format nil "rel + ldRel*~A" name)
     for abs = (format nil "abs + ldAbs*~A" name)
     for parent-name = (frame-parent f)
     do
       (if parent-name
           ;; sub frame, chain it
           (let ((par (format nil "abs + ldAbs*~A" parent-name)))
             (format stream "~&    aa_tf_qutr_~A( ~A, ~A, ~A );"
                     (if normalize "mulnorm" "mul")
                     par rel abs))
           ;; base frame, copy it
           (format stream "~&    memcpy( ~A, ~A, 7*sizeof(abs[0]) );"
                   abs rel)))
  (format stream "~&}"))

(defun emit-dot (frames &key (stream t))
  (format stream  "~&digraph {~&")
  (format stream  "~&   graph [rankdir=LR];")
  (loop for f in frames
     for name = (frame-name f)
     for parent = (frame-parent f)
     when parent
     do (format stream  "~&    ~A -> ~A;"
                parent name))
  (format stream  "~&}"))


(defun ifdef-c++ (stream string)
  (format stream "~&~%#ifdef __cplusplus~&~A~&#endif /*__cplusplus*/~%~%" string))

(defun check-frames (frames)
  ;; check indices
  (let ((hash (make-hash-table :test #'equal)))
    (loop
       for f in frames
       for i from 0
       for name =  (frame-name f)
       for parent =  (frame-parent f)
       do
         (unless name (error "Unnamed frame"))
         (setf (gethash (frame-name f) hash)
               i)
         (when parent
           (unless (gethash parent hash)
             (error "invalid parent ~A for frame ~A" parent name))))))


(defun write-frame-files (header-file source-file frames &key
                          headers
                          system-headers
                          frame-max
                          configuration-max
                          dot-file
                          normalize
                          (relative-function "rel_tf")
                          (absolute-function "abs_tf")
                          (frames-function "frames")
                          (parents-array "parents")
                          (names-array "names")
                          (axes-array "axes")
                          descriptor
                          (block-arrays t)
                          )
  (check-frames frames)
  ;; header
  (with-open-file (f header-file :direction :output :if-exists :supersede)
    (format f "~&/* AUTOGENERATED BY REFLEX FRAME GENERATOR */")
    (ifdef-c++ f "extern \"C\" {")
    (format f "~{~&#include <~A>~}" system-headers)
    (format f "~{~&#include <~A>~}" system-headers)
    (format f "~{~&#include \"~A\"~}~%" headers)
    (format f "~&extern const char *~A[];" names-array)
    (format f "~&extern const ssize_t ~A[];" parents-array)
    (format f "~&extern const double ~A[][3];" axes-array)
    (when descriptor
      (format f "~&extern const struct rfx_tf_desc ~A;" descriptor))
    (format f "~&~%/* NAME INDICES */")
    (emit-name-indices frames :stream f :max frame-max)
    (format f "~&~%/* CONFIGURATION INDICES */")
    (emit-config-indices frames :stream f :max configuration-max)
    (format f "~&~%/* Compute Relative Transforms */")
    (rel-fun-decl relative-function :stream f :block-arrays block-arrays)
    (format f ";~&")
    (format f "~&~%/* Compute Absolute Transforms */")
    (abs-fun-decl absolute-function :stream f :block-arrays block-arrays)
    (format f ";~&")
    (frame-fun-decl frames-function :stream f)
    (format f ";~&")
    (ifdef-c++ f "}"))
  ;; source
  (with-open-file (f source-file :direction :output :if-exists :supersede)
    (format f "~&/* AUTOGENERATED BY REFLEX FRAME GENERATOR */")
    (format f "~{~&#include <~A>~}" (list "amino.h" "reflex.h" header-file))
    (format f "~{~&#include \"~A\"~}" headers)
    (emit-rel-fun relative-function frames :stream f :block-arrays block-arrays)
    (emit-parents-array parents-array frames :stream f)
    (emit-names-array names-array frames :stream f)
    (emit-axes-array axes-array frames :stream f)
    (emit-abs-fun absolute-function frames :stream f :normalize normalize :block-arrays block-arrays)
    (emit-frame-fun frames-function relative-function absolute-function frames :stream f)
    (when descriptor
      (format f "~&const struct rfx_tf_desc ~A = {" descriptor)
      (format f "~&    .n_config = ~A," configuration-max)
      (format f "~&    .n_frame = ~A,"  frame-max)
      (format f "~&    .frames = ~A,"  frames-function)
      (format f "~&    .config_axes = ~A[0],"  axes-array)
      (format f "~&    .frame_name = ~A,"  names-array)
      (format f "~&    .frame_parent = ~A,"  parents-array)
      (format f "~&};"))
    )
  ;; dot
  (when dot-file
    (with-open-file (f dot-file :direction :output :if-exists :supersede)
      (emit-dot frames :stream f))))
